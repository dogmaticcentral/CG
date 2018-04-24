#!/usr/bin/python

import MySQLdb
from icgc_utils.mysql   import  *


#########################################
#########################################
def main():



	homedir = "/data/icgc"
	cancer_types = []
	for name in os.listdir(homedir):
		if os.path.isdir("/".join([homedir,name])): cancer_types.append(name)

	db     = connect_to_mysql()
	cursor = db.cursor()

	db_name =  "icgc"
	switch_to_db(cursor, db_name)
	for ct in cancer_types:
		mutations_table = ct + "_simple_somatic"
		qry = "load data local infile '%s.tsv' into table %s" % (mutations_table,mutations_table)
		search_db(cursor,qry,verbose=True)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()

