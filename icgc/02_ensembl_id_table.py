#! /usr/bin/python3

# from ensembl biomart: gene stable ID, transcript stable id

from icgc_utils.mysql import *


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

	db     = connect_to_mysql("/home/ivana/.tcga_conf")
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
