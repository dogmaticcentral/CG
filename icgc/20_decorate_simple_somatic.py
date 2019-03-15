#! /usr/bin/python3

# add info about the sequencing depth and pathogenicity
# to *simple_somatic tables for faster search

import time

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *
from config import Config

#########################################
def update_ratio_column(cursor, table):
	qry  = "update %s " % table
	qry += "set mut_to_total_read_count_ratio=mutant_allele_read_count/total_read_count "
	qry += "where total_read_count is not null and total_read_count>0"
	search_db(cursor,qry)

#########################################
def update_pathogenicity_column(cursor, table):
	qry = "select  distinct  chromosome  from %s " % table
	ret = search_db(cursor,qry)
	if not ret: return
	chroms = [line[0] for line in ret]
	for chrom in chroms:
		qry = "update  %s t, mutations_chrom_%s m " % (table, chrom)
		qry += "set t.pathogenic_estimate = m.pathogenic_estimate "
		qry += "where t.icgc_mutation_id=m.icgc_mutation_id "
		search_db(cursor,qry)

########################################
def decorate(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:

		print("====================")
		print("decorating ", table, os.getpid())

		time0 = time.time()

		###############
		print("\t\t updating ratio column in {} ".format(table))
		update_ratio_column(cursor, table)

		###############
		print("\t\t updating pathogenicity column in {} ".format(table))
		update_pathogenicity_column(cursor, table)

		time1 = time.time()
		print("\t\t %s done in %.3f mins" % (table, float(time1-time0)/60)), os.getpid()

	cursor.close()
	db.close()

	return

#########################################
#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in search_db(cursor,qry)]
	cursor.close()
	db.close()

	number_of_chunks = 20
	parallelize (number_of_chunks, decorate, tables, [])


	return



#########################################
if __name__ == '__main__':
	main()
