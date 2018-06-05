#! /usr/bin/python
import subprocess


from icgc_utils.common_queries import *

def protein_mutations (cursor, gene_symbol, p53_status='wt'):

	# a different query - I want this over all tables
	qry =  "select g1.gene_symbol,  g2.gene_symbol, s1.icgc_donor_id, s1.submitted_sample_id  "
	qry += "from mutation2gene g1, mutation2gene g2,  %s s1,  %s s2  " % (somatic_table, somatic_table)
	qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
	qry += "and s1.icgc_mutation_id=g1.icgc_mutation_id and g1.gene_symbol='%s' " % gene_symbol
	qry += "and s2.icgc_mutation_id=g2.icgc_mutation_id and g2.gene_symbol='%s' " % 'TP53'
	qry += "and s1.pathogenic_estimate=1 and s1.reliability_estimate=1 "
	qry += "and s2.pathogenic_estimate=1 and s2.reliability_estimate=1 "

	return [r[0] for r in search_db(cursor,qry)]


def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	gene = 'RPL5'

	for p53_status in ['wt', 'mutated']:
		mutations = protein_mutations (cursor, gene, p53_status=p53_status)
	
		# clustering input

		# clustering run


		# parse clustering output

		# pymol input



	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
