#! /usr/bin/python

from icgc_utils.common_queries import *

verbose = False

def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")
	total_donors = 0
	total_rpl5 = 0
	total_tp53 = 0
	total_cooc = 0
	for table in tables:
		patients_with_muts_in_gene = patients_per_gene_breakdown(cursor, table)
		if patients_with_muts_in_gene.get('RPL5', 0)==0: continue

		print "================================="
		print table
		donors = len(get_donors(cursor, table))
		print "donors: ", donors
		total_donors += donors
		for gene in ['RPL5', 'TP53']:
			print gene, patients_with_muts_in_gene.get(gene, 0)
		total_rpl5 += patients_with_muts_in_gene.get('RPL5',0)
		total_tp53 += patients_with_muts_in_gene.get('TP53',0)
		#co_ocurrence = 	co_ocurrence_raw(cursor, table, 'RPL5', 'TP53')
		#if not co_ocurrence:
		#	print "no co-ocurrence"
		#else:
		#	for line in co_ocurrence:
		#		print line
		#exit()
		cooc = co_ocurrence_count(cursor, table, 'RPL5', 'TP53')
		print "co-ocurrence:", cooc
		total_cooc += cooc

		# use littler to evaluate


	print
	print "================================="
	print "total donors:", total_donors
	print "        rpl5:", total_rpl5
	print "        tp53:", total_tp53
	print "        cooc:", total_cooc
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
