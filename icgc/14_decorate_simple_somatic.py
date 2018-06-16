#! /usr/bin/python

# add info about the sequencing depth and pathogenicity
# to *simple_somatic tables for faster search

import time

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *

#########################################
def add_columns(cursor, table):

	# mut_to_total_read_count_ratio
	column_name = "mut_to_total_read_count_ratio"
	if not column_exists(cursor, "icgc", table, column_name):
		qry = "ALTER TABLE %s ADD  %s float default 0.0" % (table, column_name)
		search_db(cursor,qry, verbose=True)
	for column_name in  ["pathogenic_estimate","is_missense"]:
		if not column_exists(cursor, "icgc", table, column_name):
			qry = "ALTER TABLE %s ADD  %s boolean  default 0" % (table, column_name)
			search_db(cursor,qry, verbose=True)



#########################################
def update_ratio_column(cursor, table):
	qry  = "update %s " % table
	qry += "set mut_to_total_read_count_ratio=mutant_allele_read_count/total_read_count "
	qry += "where total_read_count is not null and total_read_count>0"
	search_db(cursor,qry, verbose=True)

#########################################
def update_pathogenicity_column(cursor, table):
	qry = "select  distinct  chromosome  from %s " % table
	chroms =[ret[0] for ret in  search_db(cursor,qry)]
	for chrom in chroms:
		qry = "update  %s t, mutations_chrom_%s m " % (table, chrom)
		qry += "set t.pathogenic_estimate = m.pathogenic_estimate "
		qry += "where t.icgc_mutation_id=m.icgc_mutation_id "
		search_db(cursor,qry)

#########################################
def update_missense_column(cursor, table):
	qry = "select  distinct  chromosome  from %s " % table
	chroms =[ret[0] for ret in  search_db(cursor,qry)]
	for chrom in chroms:

		qry = "update  %s t, mutations_chrom_%s m " % (table, chrom)
		qry += "set t.is_missense = 1 "
		qry += "where m.consequence like '%missense%' "
		qry += "where t.icgc_mutation_id=m.icgc_mutation_id "
		search_db(cursor,qry)



########################################
def decorate(tables, other_args):

	db     = connect_to_mysql()
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:

		time0 = time.time()

		print "===================="
		print "decorating ", table, os.getpid()

		add_columns(cursor, table)
		#update_ratio_column(cursor, table)

		###############
		##print "\t\t updating pathogenicity column in %s " % table
		##update_pathogenicity_column(cursor, table)
		##time1 = time.time()
		##print ("\t\t %s done in %.3f mins" % (table, float(time1-time0)/60)), os.getpid()

		###############
		print "\t\t updating missense column in %s " % table
		qry = "ALTER TABLE %s drop column is_missense " % table
		search_db(cursor,qry,verbose=True)
		#update_missense_column(cursor, table)
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
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	number_of_chunks = 8
	parallelize (number_of_chunks, decorate, tables, [])


	return



#########################################
if __name__ == '__main__':
	main()
