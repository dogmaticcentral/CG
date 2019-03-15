#! /usr/bin/python3

from icgc_utils.common_queries import *
from icgc_utils.processes import *
from config import Config


#########################################
def add_columns(cursor, table):

	# mut_to_total_read_count_ratio
	column_name = "reliability_estimate"
	if not column_exists(cursor, "icgc", table, column_name):
		qry = "ALTER TABLE %s ADD  %s boolean  default 0" % (table, column_name)
		search_db(cursor,qry, verbose=True)


#########################################
def add_reliability_annotation(tables, other_args):
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:
		print(table, "running")
		qry  = "update %s set reliability_estimate = 1 " % table
		# for some data sets we do not have the read count info
		qry += "where total_read_count is null or "
		qry += "(mutant_allele_read_count>=10 and mut_to_total_read_count_ratio>=0.2)"
		search_db(cursor,qry)
	cursor.close()
	db.close()


#########################################
#########################################
def main():

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	for table in tables:
		print("checking/adding reliability column to", table)
		add_columns(cursor, table)
	cursor.close()
	db.close()

	print("adding annotation")

	number_of_chunks = 20
	parallelize(number_of_chunks, add_reliability_annotation, tables, [])


	return



#########################################
if __name__ == '__main__':
	main()
