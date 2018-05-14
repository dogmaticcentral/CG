#! /usr/bin/python
import time

from icgc_utils.mysql   import  *
from icgc_utils.processes import *
from random import shuffle

def make_map_table(cursor, db_name, table_name):

	if check_table_exists(cursor, db_name, table_name): return

	switch_to_db (cursor, db_name)

	qry = ""
	qry += "  CREATE TABLE  %s (" % table_name
	qry += "     id INT not null AUTO_INCREMENT, "
	qry += "  	 icgc_mutation_id VARCHAR (20) NOT NULL, "
	qry += "	 gene_symbol VARCHAR (30) NOT NULL, "
	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print qry
	print rows

#########################################
def store (cursor, mut_id, gene_symbols):
	for symbol in gene_symbols:
		qry  = "insert into mutation2gene (icgc_mutation_id,gene_symbol) "
		qry += "values ('%s','%s')" % (mut_id, symbol)
		if search_db(cursor,qry):
			search_db(cursor,qry,verbose=True)

	return

#########################################
def ens2hgnc(cursor,ensids):
	symbols = set([])
	for ensid in ensids:
		qry  = "select approved_symbol from hgnc "
		qry += "where locus_group='protein-coding gene' "
		qry += "and (ensembl_gene_id='%s' or ensembl_gene_id_by_hgnc='%s')" % (ensid,ensid)
		ret = search_db(cursor, qry)
		if not ret: continue
		if len(ret)>1:
			print "nonunique symbol for %s (?)" % ensid,
			print ",".join([r[0] for r in ret])
			#continue
		symbols.add(ret[0][0])

	return symbols

#########################################
def transcr2gene(cursor, transcrids):

	switch_to_db(cursor,"homo_sapiens_core_91_38")

	geneids = set([])
	for trscid in transcrids:

		qry  = "select g.stable_id from gene g, transcript t "
		qry += "where g.gene_id = t.gene_id and t.stable_id='%s'" % trscid
		ret = search_db(cursor,qry)
		if not ret: continue
		geneids.add(ret[0][0])

	switch_to_db(cursor,"icgc")
	return geneids

#########################################
def store_maps(chromosomes, other_args ):
	db     = connect_to_mysql()
	cursor = db.cursor()
	switch_to_db(cursor, "icgc")
	for chrom in chromosomes:
		time0 = time.time()
		print "===================="
		print "maps for ", chrom, "pid:", os.getpid()

		qry  = "select m.icgc_mutation_id, l.gene_relative, l.transcript_relative "
		qry += "from mutations_chrom_%s m, locations_chrom_%s l " % (chrom, chrom)
		qry += "where m.start_position=l.position and (l.gene_relative is not null or l.transcript_relative is not null)"
		ret = search_db(cursor, qry, verbose=True)
		if not ret:
			print "(?) no ret for "
			print  qry
			exit()
		for mut_id, gene, transcr in ret:
			gene_found = False
			if gene and gene != "":
				geneids = set ([])
				for geneloc in gene.split(";"):
					ensid, loc = geneloc.split(":")
					if loc == "intragenic": geneids.add(ensid)

				if len(geneids)>0:
					#print mut_id, ensids
					gene_found = True
					store (cursor, mut_id, ens2hgnc(cursor,geneids))
			if not gene_found and transcr and transcr != "":
				transcrids = set ([])
				for transcrloc in transcr.split(";"):
					transcrid, loc = transcrloc.split(":")
					transcrids.add(transcrid)
					geneids = transcr2gene(cursor, transcrids)
					store (cursor, mut_id, ens2hgnc(cursor,geneids))
		time1 = time.time()
		print "chrom ", chrom, "done in %.3f mins" % (float(time1-time0)/60)
	cursor.close()
	db.close()

	return


#########################################
#########################################
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()
	#########################
	make_map_table(cursor, "icgc", "mutation2gene")

	chromosomes = [str(i) for i in range(1,23)] + ["X","Y"]
	shuffle(chromosomes)
	number_of_chunks = 8
	parallelize (number_of_chunks, store_maps, chromosomes, [])

	return

#########################################
if __name__ == '__main__':
	main()
