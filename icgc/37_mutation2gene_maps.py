#! /usr/bin/python3
#
# This source code is part of icgc, an ICGC processing pipeline.
# 
# Icgc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Icgc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#
import time

from icgc_utils.mysql import *
from icgc_utils.processes import *
from random import shuffle
from config import Config

#########################################
def make_map_table(cursor, db_name, table_name):

	if check_table_exists(cursor, db_name, table_name): return

	switch_to_db (cursor, db_name)
	# we cannot use icgc_mutation_id bcs the same location can correspond to multiple genes
	qry = ""
	qry += "  CREATE TABLE  %s (" % table_name
	qry += "     id INT not null AUTO_INCREMENT, "
	qry += "  	 icgc_mutation_id VARCHAR (20) NOT NULL, "
	qry += "	 gene_symbol VARCHAR (30) NOT NULL, "
	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	search_db(cursor, qry)

	# searching for priors will take the longest time though
	qry = "create index mut_idx on mutation2gene (icgc_mutation_id)"
	search_db(cursor, qry, verbose=True)


#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./icgc_utils/kernprof.py -l 37_....py
# followed by
# python3 -m line_profiler 37_....py.lprof
# see here https://github.com/rkern/line_profiler#line-profiler
# the reason I am using local kernprof.py is that I don't know where pip
# installed its version (if anywhere)
# @profile
#########################################
def store (cursor, mut_id, gene_symbols):

	# check existing
	qry = "select gene_symbol from mutation2gene "
	qry += "where icgc_mutation_id='%s'" % mut_id
	ret = search_db(cursor,qry)
	if not ret:
		existing = None
	else:
		existing = set([r[0] for r in ret])
	if existing:
		gene_symbols = list(set(gene_symbols).difference(existing))
	#insert new
	for symbol in gene_symbols:
		qry  = "insert into mutation2gene (icgc_mutation_id,gene_symbol) "
		qry += "values ('%s','%s')" % (mut_id, symbol)
		if search_db(cursor,qry):
			search_db(cursor,qry,verbose=True)
			exit()
	return

#########################################
def ens2hgnc(cursor,ensids):
	symbols = set([])
	for ensid in ensids:
		# note that we are filtering for protein-coding genes here:
		qry  = "select approved_symbol from hgnc "
		qry += "where locus_group='protein-coding gene' "
		qry += "and (ensembl_gene_id='%s' or ensembl_gene_id_by_hgnc='%s')" % (ensid,ensid)
		ret = search_db(cursor, qry)
		if not ret: continue
		if len(ret)>1:
			print("nonunique symbol for %s (?)" % ensid)
			print(",".join([r[0] for r in ret]))
			print(qry())
			exit()
		symbols.add(ret[0][0])

	return symbols


#########################################
def transcr2gene(cursor, transcrids):

	geneids = []
	qry  = "select gene from ensembl_ids "
	qry += "where transcript in (%s)" % (",".join(["'%s'"%tr for tr in transcrids]))
	ret = search_db(cursor,qry)
	if ret: geneids = [r[0] for r in ret]
	return geneids


#########################################
def report_progress(chrom, ct, no_rows, time0):
	if (ct%10000>0): return
	print("%30s   %6d lines out of %6d  (%d%%)  %d min" % \
	(chrom, ct, no_rows, float(ct)/no_rows*100, float(time.time()-time0)/60))
	return

#########################################
def store_maps(chromosomes, other_args ):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, "icgc")
	for chrom in chromosomes:
		time0 = time.time()
		print("====================")
		print("maps for ", chrom, "pid:", os.getpid())

		qry = "select m.icgc_mutation_id, l.gene_relative, l.transcript_relative "
		qry += "from mutations_chrom_%s m, locations_chrom_%s l " % (chrom, chrom)
		qry += "where m.start_position=l.position and (l.gene_relative is not null or l.transcript_relative is not null)"
		ret = search_db(cursor, qry, verbose=True)

		if not ret:
			print("(?) no ret for ")
			print(qry)
			exit()
		no_rows = len(ret)
		print("chrom ", chrom, "number of rows", no_rows)

		ct = 0
		for mut_id, gene, transcr in ret:
			ct += 1
			report_progress(chrom, ct, no_rows, time0)
			gene_found = False
			if gene and gene != "":
				geneids = set ([])
				for geneloc in gene.split(";"):
					ensid, loc = geneloc.split(":")
					if loc == "intragenic": geneids.add(ensid)
				if len(geneids)>0:
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
		print("chrom ", chrom, "done in %.3f mins" % (float(time1-time0)/60))
	cursor.close()
	db.close()

	return


#########################################
#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	make_map_table(cursor, "icgc", "mutation2gene")

	chromosomes = [str(i) for i in range(1,23)] + ["X","Y"]
	shuffle(chromosomes)

	number_of_chunks = 10
	parallelize (number_of_chunks, store_maps, chromosomes, [])

	return

#########################################
if __name__ == '__main__':
	main()
