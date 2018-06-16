#! /usr/bin/python

# add info about the reliability to mutations_chrom_* tables
# for faster search

import time

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *
from random import shuffle

#########################################
def add_columns(cursor, table):

	# mut_to_total_read_count_ratio
	column_name = "reliability_estimate"
	if not column_exists(cursor, "icgc", table, column_name):
		qry = "ALTER TABLE %s ADD  %s boolean  default 0" % (table, column_name)
		search_db(cursor,qry, verbose=True)


#########################################
# this will be set to 1 if we find support for this variant anywhere
def update_reliability_column_in_mutation_tables(cursor, table):
	qry = "select  distinct  chromosome  from %s " % table
	chroms =[ret[0] for ret in  search_db(cursor,qry)]
	for chrom in chroms:
		qry = "update  %s t, mutations_chrom_%s m " % (table, chrom)
		qry += "set m.reliability_estimate = 1 "
		qry += "where t.icgc_mutation_id=m.icgc_mutation_id "
		qry += "and  (t.total_read_count is null or "
		# since the depths for the matched samples are not given, we are using these quite strict rules ...
		qry += "(t.mutant_allele_read_count>3 and t.mut_to_total_read_count_ratio>0.2) )"
		search_db(cursor,qry)

########################################
def decorate_mutations(tables, other_args):

	db     = connect_to_mysql()
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:

		time0 = time.time()
		print "===================="
		print "using annotations from  ", table,  "to add reliability info to  mutations; pid:", os.getpid()
		update_reliability_column_in_mutation_tables(cursor, table)
		time1 = time.time()
		print ("\t\t %s done in %.3f mins" % (table, float(time1-time0)/60)), os.getpid()


	cursor.close()
	db.close()

	return

#########################################
#########################################
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()
	chromosomes = [str(i) for i in range(1,23)] + ["X","Y"]
	switch_to_db(cursor,"icgc")
	for chromosome in chromosomes:
		table = "mutations_chrom_%s"%chromosome
		print "===================="
		print "checking/adding reliability column to", table
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
