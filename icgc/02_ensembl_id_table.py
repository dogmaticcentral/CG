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

# from ensembl biomart: gene stable ID, transcript stable id

from icgc_utils.mysql import *
from config import Config

#########################################
# make one-on-one ENST to ENSG translation table
def make_ensembl_ids_table(cursor, db_name, ens_ids_table):
	if check_table_exists (cursor, db_name, ens_ids_table): return
	switch_to_db (cursor, db_name)
	qry = ""
	qry += "  CREATE TABLE  %s (" % ens_ids_table
	charlen = 20
	for name in ['transcript', 'gene', 'canonical_transcript']:
		qry += " %s VARCHAR(%d) NOT NULL ," % (name, charlen)

	qry += "	 PRIMARY KEY (transcript) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)
	return

#########################################
#########################################
def main():

	ens_id_file = "/storage/databases/ensembl-94/ensembl_gene2trans_stable.tsv"

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	db_name = "icgc"

	make_ensembl_ids_table(cursor, db_name, "ensembl_ids")

	qry = "load data local infile '%s' into table %s" % (ens_id_file,"ensembl_ids")
	search_db(cursor, qry, verbose=True)


	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
