#! /usr/bin/python3

# add info about the reliability to mutations_chrom_* tables
# for faster search

import time

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *
from config import Config

#########################################
def add_columns(cursor, table):

	# mut_to_total_read_count_ratio
	column_name = "reliability_estimate"
	if not column_exists(cursor, "icgc", table, column_name):
		qry = "ALTER TABLE %s ADD  %s boolean  default 0" % (table, column_name)
		search_db(cursor,qry, verbose=True)


#########################################
def update_reliability_column_in_mutation_tables(cursor, table):
	chromosomes = [str(i) for i in range(1,23)] + ["X","Y"]
	for c in chromosomes:
		# is there _any_ reliable source for this mutation?
		# if so, set reliability to 1
		qry  = "update  %s s, mutations_chrom_%s m " % (table, c)
		qry += "set m.reliability_estimate=1  "
		qry += "where m.reliability_estimate=0 and s.reliability_estimate=1 "
		qry += "and s.icgc_mutation_id = m.icgc_mutation_id "
		search_db(cursor,qry, verbose=False)


########################################
def decorate_mutations(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:

		time0 = time.time()
		print("====================")
		print("using annotations from  ", table, "to add reliability info to  mutations; pid:", os.getpid())
		update_reliability_column_in_mutation_tables(cursor, table)
		time1 = time.time()
		print(("\t\t %s done in %.3f mins" % (table, float(time1-time0)/60)), os.getpid())

	cursor.close()
	db.close()

	return

#########################################
#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	chromosomes = [str(i) for i in range(1,23)] + ["X","Y"]
	switch_to_db(cursor,"icgc")
	for chromosome in chromosomes:
		table = "mutations_chrom_%s"%chromosome
		print("====================")
		print("checking/adding reliability column to", table)
		add_columns(cursor, table)

	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	number_of_chunks = 8
	parallelize(number_of_chunks, decorate_mutations, tables, [])


	return



#########################################
if __name__ == '__main__':
	main()
