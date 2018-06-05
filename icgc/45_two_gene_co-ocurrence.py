#! /usr/bin/python
import subprocess


from icgc_utils.common_queries import *
from scipy import stats
from math import exp, log

verbose = False
def myfisher(donors, gene_1_mutated, other_mutated, cooc):
	# https://en.wikipedia.org/wiki/Hypergeometric_distribution
	# https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html
	# In probability theory and statistics, the cumulative distribution function
	# (CDF, also cumulative density function) of a real-valued random variable X, or just distribution function of X,
	#  evaluated at x,
	# is the probability that X will take a value less than or equal to x
	p_smaller = stats.hypergeom.cdf(cooc, donors, gene_1_mutated, other_mutated)
	if cooc>0: # the less or equal in both cases is the reason why they do not add up to 1
		p_greater = 1 -  stats.hypergeom.cdf(cooc-1, donors, gene_1_mutated, other_mutated)
	else:
		p_greater = 1.0
	return p_smaller, p_greater


#####################################
def fisher(donors, gene_1_mutated, other_mutated, cooc):
	a = donors - gene_1_mutated - other_mutated + cooc # both wt (we subtracted the overlap twice)
	b = gene_1_mutated - cooc  # p53  mutated adn rpl5 wt
	c = other_mutated - cooc # rpl5 mutated and p53 wt
	d = cooc                 # both mutated
	[odds, pval_lt] = stats.fisher_exact([[a, b], [c, d]], "less")
	[odds, pval_gt] = stats.fisher_exact([[a, b], [c, d]], "greater")


	return pval_lt, pval_gt



###################################
# ## file:///home/ivana/Dropbox/Sinisa/ribosomal/html/the_curious_case_of_rpl22.html
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	gene_1 = 'TP53'
	#gene_1 = 'BRCA1'
	other_genes = ['RPL5']
	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]


	#########################
	switch_to_db(cursor,"icgc")
	total_donors = 0
	total_other = 0
	total_gene_1 = 0
	total_cooc = 0
	#tables = ['UCEC_simple_somatic']

	for table in tables:
		tumor_short = table.split("_")[0]
		patients_with_muts_in_gene = patients_per_gene_breakdown(cursor, table)
		if patients_with_muts_in_gene.get(gene_1,0)==0: continue
		no_mutant = True
		for gene_2 in other_genes:
			if patients_with_muts_in_gene.get(gene_2,0)==0: continue
			no_mutant = False
			break
		if no_mutant: continue

		print "================================="
		print table
		donors = len(get_donors(cursor, table))
		print "donors: ", donors
		total_donors += donors
		for gene in [gene_1]+other_genes:
			print gene, patients_with_muts_in_gene.get(gene, 0)

		gene_1_mutated  = patients_with_muts_in_gene.get(gene_1,0)
		other_mutated = patients_with_muts_in_gene_group(cursor, table, other_genes)
		total_gene_1   += gene_1_mutated
		total_other  += other_mutated

		cooc = co_ocurrence_w_group_count(cursor, table, gene_1, other_genes)

		p_smaller, p_bigger = myfisher(donors, gene_1_mutated, other_mutated, cooc)


		#pval_lt, pval_gt = fisher(donors, gene_1_mutated, other_mutated, cooc)

		expected = float(gene_1_mutated)/donors*other_mutated
		print "co-ocurrence:", cooc
		print "    expected: %.1f" % expected
		print "   p_smaller: %.2f" % p_smaller
		print "    p_bigger: %.2f" % p_bigger
		print
		total_cooc += cooc



	#pval_lt, pval_gt = fisher(total_donors, total_gene_1, total_other, total_cooc)
	p_smaller, p_bigger = myfisher(total_donors, total_gene_1, total_other, total_cooc)
	print
	print "================================="
	print other_genes
	print "total donors:", total_donors
	print "        other:", total_other
	print "          %s: %d" % (gene_1, total_gene_1)
	print "        cooc:", total_cooc
	print "    expected: %.1f" % (float(total_gene_1)/total_donors*total_other)
	#print "     pval_lt: %.1e" % pval_lt
	#print "     pval_gt: %.1e" % pval_gt
	print "   p_smaller: %.1e" % p_smaller
	print "    p_bigger: %.1e" % p_bigger
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
