#! /usr/bin/python


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
def copy_reliability_annotation(tables, other_args):
	db     = connect_to_mysql()
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	chromosomes = [str(i) for i in range(1,23)] + ["X","Y"]
	for table in tables:
		print table, "running"
		for c in chromosomes:
			qry  = "update  %s s, mutations_chrom_%s m " % (table, c)
			qry += "set s.reliability_estimate = m.reliability_estimate "
			qry += "where s.icgc_mutation_id = m.icgc_mutation_id"
			search_db(cursor,qry, verbose=False)
		print table, "done"
	cursor.close()
	db.close()

#########################################
#########################################
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	for table in tables:
		print "===================="
		print "checking/adding reliability column to", table
		add_columns(cursor, table)

	cursor.close()
	db.close()

	print "copying annotation"

	number_of_chunks = 8
	parallelize(number_of_chunks, copy_reliability_annotation, tables, [])

	return



#########################################
if __name__ == '__main__':
	main()
