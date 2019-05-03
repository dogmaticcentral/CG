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

# The UCSC Genome Browser database: 2019 update.: pubmed id 30407534
# To keep in mind:  all start coordinates in UCSC database are 0-based, not 1-based:
# http://genome.ucsc.edu/FAQ/FAQtracks#tracks1
# it looks like -1 frame means it is UTR (https://www.biostars.org/p/160487/)

from config import Config
from icgc_utils.mysql import *

#########################################
#  "select g.name2, g.name, g.strand, g.txStart, g.txEnd, g.cdsStart, g.cdsEnd, g.exonCount, g.exonStarts, g.exonEnds "
def make_gene_coords_table(cursor, db_name, coords_table):
	switch_to_db (cursor, db_name)

	if check_table_exists (cursor, db_name, coords_table):
		search_db(cursor, "drop table " + coords_table)

	qry = ""
	qry += "  CREATE TABLE  %s (" % coords_table
	qry += "	 ens_transcript_id  VARCHAR (20) NOT NULL, "
	qry += "  	 ens_gene_id VARCHAR (20) NOT NULL, "
	qry += "	 strand  CHAR (1) NOT NULL, "

	qry += "	 tx_start INT  NOT NULL, "
	qry += "	 tx_end  INT NOT NULL, "

	qry += "	 cds_start INT  NOT NULL, "
	qry += "	 cds_end  INT NOT NULL, "

	qry += "	 exon_count  INT NOT NULL, "

	qry += "     exon_starts  TEXT, "
	qry += "     exon_ends TEXT, "
	qry += "     exon_frames TEXT, "

	qry += "	 PRIMARY KEY (ens_transcript_id) "
	qry += ") ENGINE=MyISAM"

	search_db(cursor, qry)
	return



#########################################
def main():

	chromosomes = [str(x) for x in range(1,23)] + ["X", "Y"]
	assembly = "hg19"

	for dependency in [Config.mysql_conf_file, Config.ucsc_mysql_conf_file]:
		if not os.path.exists(dependency):
			print(dependency, "not found")
			exit()

	local_db = connect_to_mysql(Config.mysql_conf_file)
	local_cursor = local_db.cursor()
	switch_to_db(local_cursor,'icgc')

	ucsc_db     = connect_to_mysql(Config.ucsc_mysql_conf_file)
	ucsc_cursor = ucsc_db.cursor()
	switch_to_db(ucsc_cursor, assembly)

	home = os.getcwd()
	# make a workdir and move there
	workpath = "{}/tsvs/ucsc".format(home)
	if not os.path.exists(workpath): os.makedirs(workpath)
	os.chdir(workpath)

	for chrom in chromosomes:


		fnm  = "chrom_%s.tsv" % chrom
		outf = open(fnm,"w")
		print("downloading data for", assembly, chrom)
		# note we will take only genes that have a known associated peptide
		# tx stands for transcription; ensGene schema:
		# http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?hgsid=6533228_Aczpyj65LpzrrJR9FAo04iooB2aa&hgta_doSchemaDb=hg19&hgta_doSchemaTable=ensGene
		# name is ENST, name 2 is ENSG
		qry  = "select g.name, g.name2, g.strand, g.txStart, g.txEnd, g.cdsStart, g.cdsEnd, "
		qry += "g.exonCount, g.exonStarts, g.exonEnds, g.exonFrames "
		qry += "from ensGene  g right join ensPep p on g.name=p.name "
		qry += "where g.chrom='chr%s'" % chrom
		rows = search_db(ucsc_cursor,qry)

		print("storing to tsv " + fnm)
		for row in rows:
			#[gene_name, strand, rfrom, rto] = row
			outf.write("\t".join([str(r) for r in row]) + "\n")
		outf.close()

		local_table = "coords_chrom_%s"%chrom
		print("making gene coords table for " + chrom)
		make_gene_coords_table(local_cursor, "icgc", local_table)

		print("loading gene coords table for " + chrom)
		qry = "load data local infile '%s/%s' into table %s" % (workpath, fnm, local_table)
		search_db(local_cursor,qry)

		print("indexing " + local_table)
		qry  = "create index tx_range_idx on %s (tx_start,tx_end)" % local_table
		search_db(local_cursor,qry)


	ucsc_cursor.close()
	ucsc_db.close()


	local_cursor.close()
	local_db.close()
	return True


#########################################
if __name__ == '__main__':
	main()


