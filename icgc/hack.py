#! /usr/bin/python
import subprocess
import time, re

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *

#########################################
def main():

	# divide by cancer types, because I have duplicates within each cancer type
	# that I'll resolve as I go, but I do not want the threads competing)
	db     = connect_to_mysql()
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in search_db(cursor,qry)]

	for table in tables:
		if column_exists (cursor, 'icgc', table,'reliability_estimate'): continue
		print table
		qry  = "alter table icgc.%s " % table
		#qry += "modify column control_genotype varchar(430) default null"
		qry  += "add  reliability_estimate  boolean  default 0"
		search_db(cursor,qry, verbose=True)

		#search_db(cursor,qry,verbose=True)
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
