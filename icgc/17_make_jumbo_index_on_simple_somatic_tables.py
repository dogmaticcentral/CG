#!/usr/bin/python

import MySQLdb

from config import Config
from icgc_utils.mysql   import  *

#########################################
#########################################
def main():


	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")

	# indices on simple somatic temp
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	for table in tables:
		print(table)
		qry   = "create index mega_idx on %s " % table
		qry += " (icgc_mutation_id, icgc_donor_id, icgc_specimen_id, icgc_sample_id)"
		search_db(cursor,qry,verbose=True)
		# while at it, let's add this one too (to be used in 18_cleanup_duplicate_specimens)
		qry   = "create index donor_idx on %s (icgc_donor_id)" % table
		search_db(cursor,qry,verbose=True)


	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()

