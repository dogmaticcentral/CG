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
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in search_db(cursor,qry)]

	switch_to_db(cursor,'icgc')

	for table in tables:
		print table
		qry  = "select icgc_donor_id from icgc.%s " % table
		qry += "where icgc_mutation_id ='MU1800992'"
		search_db(cursor,qry)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
