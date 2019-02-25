#!/usr/bin/python3

import MySQLdb

from config import Config
from icgc_utils.mysql   import  *

#########################################
#########################################
def main():

	homedir = Config.data_home_local
	cancer_types = []
	for name in os.listdir(homedir):
		if os.path.isdir("/".join([homedir,name])): cancer_types.append(name)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")
	cons_vocab = set([])
	for ct in cancer_types:
		mutations_table = ct + "_simple_somatic_temp"

		qry = "select distinct consequence_type from %s" % mutations_table
		csq = [ret[0] for ret in search_db(cursor,qry)]
		cons_vocab = cons_vocab.union(csq)

	for kwd in cons_vocab:
		print(kwd)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()



