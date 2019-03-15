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
	qry += "where table_schema='icgc' and table_name like '%_specimen'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	for table in tables:
		print(table)
		qry   = "alter table %s " % table
		qry += " modify column icgc_specimen_id varchar(50) "
		search_db(cursor,qry,verbose=True)


	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()

