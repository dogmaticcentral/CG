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

# from HUGO gene nomenclature committee
# https://www.genenames.org/download/custom/

# header names:
# head -n1 /data/hgnc/hgnc_name_res.tsv \
# | sed 's/\t/\n/g' | awk 'ct +=1 {printf "%d ", ct; print}'

# locus group:
# | protein-coding gene |
# | non-coding RNA      |
# | pseudogene          |
# | other               |
# | phenotype           |

from config import Config
from icgc_utils.mysql import *

header_names = ["hgnc_id", "approved_symbol", "approved_name",
				"locus_group", "synonyms", "chromosome",
				"ensembl_gene_id_by_hgnc", "refseq_ids", "uniprot_ids",
				"ensembl_gene_id"]

#########################################
def make_hgnc_table(cursor, db_name, hgnc_table):

	if check_table_exists (cursor, db_name, hgnc_table): return
	switch_to_db (cursor, db_name)

	qry = ""
	qry += "  CREATE TABLE  %s (" % hgnc_table
	qry += "     id INT NOT NULL, "
	for name in header_names:
		if name in ['approved_name', 'synonyms', 'refseq_ids']:
			charlen = 150
		elif name == 'uniprot_ids':
			charlen = 300
		elif name in ['approved_symbol', 'chromosome']: # chrom can have annotation such as "not on reference assembly"
			charlen = 30
		else:
			charlen = 20
		qry += " %s VARCHAR(%d)," % (name, charlen)

	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)
	return

#########################################
def strip_arm_annotation(chrom_address):
	if "p" in chrom_address:
		return chrom_address.split("p")[0]
	if "q" in chrom_address:
		return chrom_address.split("q")[0]
	return chrom_address

#########################################
#########################################
def main():

	hgncfile = "/storage/databases/hgnc/hgnc_name_res.tsv"

	tmp_outfile = "hgnctmp.tsv"
	ct = 0
	outf = open (tmp_outfile, "w")
	with open(hgncfile, "r") as inf:
		headers = inf.readline().rstrip("\n").split("\t")
		chromosome_column = headers.index('Chromosome')
		print(chromosome_column)
		for line in inf:
			fields = line.rstrip("\n").split("\t")
			fields[chromosome_column] = strip_arm_annotation(fields[chromosome_column])
			ct += 1
			outfields = [str(ct)] + fields
			outf.write("\t".join(outfields)+"\n")
	outf.close()

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	db_name = "icgc"

	make_hgnc_table(cursor, db_name, "hgnc")
	qry = "load data local infile '%s' into table %s" % ("hgnctmp.tsv","hgnc")
	search_db(cursor,qry,verbose=True)
	os.remove(tmp_outfile)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
