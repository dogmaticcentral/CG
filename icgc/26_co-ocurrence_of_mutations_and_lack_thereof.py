#! /usr/bin/python

from icgc_utils.common_queries import *

verbose = False

def co_ocurrence(cursor,table,gene):

	total_rpl5 += patients_with_muts_in_gene.get('RPL5',0)
	total_tp53 += patients_with_muts_in_gene.get('TP53',0)
	cooc = co_ocurrence_count(cursor, table, 'RPL5', 'TP53')
	print "co-ocurrence:", cooc
	total_cooc += cooc

# use littler to evaluate




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

	# the total sums are per gene, over tumors that have a mutation in this gene
	total_donors = {}
	total = {}
	total_tp53 = {} # total tp53 mutated in tumors that have gene of interest mutated
	total_cooc = {}
	for line in search_db(cursor,"select distinct(gene_symbol) from mutation2gene"):
		gene = line[0]
		total_donors[gene] = 0
		total[gene] = 0
		total_tp53[gene] = 0
		total_cooc[gene] = 0

	for table in tables:
		print "================================="
		print table
		patients_with_muts_in_gene = patients_per_gene_breakdown(cursor, table)
		donors = len(get_donors(cursor, table))
		print "donors: ", donors
		for gene, number_of_patients in patients_with_muts_in_gene.iteritems():
			if gene=='TP53': continue
			total_donors[gene] += donors
			total_tp53[gene]   += patients_with_muts_in_gene.get('TP53',0)
			total[gene] += number_of_patients
			total_cooc[gene] += co_ocurrence_count(cursor, table, gene, 'TP53')

	outf = open("coocurrence.table","w")
	for gene, donors in total_donors.iteritems():
		if donors<10: continue
		outf.write(" %10s  %5d  %4d  %4d  %4d \n" % (gene, donors, total_tp53[gene], total[gene], total_cooc[gene]))
	outf.close()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
