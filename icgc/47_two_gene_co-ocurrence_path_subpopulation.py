#! /usr/bin/python

# limiting ourselves to pathogenic population in gene_2
# - do missense mutations have a greater tendency to co-occur with wt tp53?

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
	gene_2  = 'RPL11'
	chromosome = find_chromosome(cursor, gene_2)
	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]


	#########################
	switch_to_db(cursor,"icgc")
	total_donors = 0
	total_p53path = 0
	total_miss = 0
	total_cooc = 0

	write_to_file =  False

	if write_to_file:
		outf = open("{}_{}_cooccurrence.tsv".format(gene_1, gene_2),"w")
		outf.write("\t".join(['cancer','donors', "mutations in %s"%gene_1,
			                "mutations in %s"%gene_2, 'cooccurrence','expected',
			                'p smaller', 'p bigger'])+"\n")

	for table in tables:
		print "================================="
		print table

		tumor_short = table.split("_")[0]

		qry  = "select s.icgc_donor_id,  s.icgc_mutation_id, s.icgc_specimen_id "
		qry += "from mutation2gene g,  %s s " % table
		qry += "where g.gene_symbol='%s' " % gene_2
		qry += "and g.icgc_mutation_id = s.icgc_mutation_id "
		qry += "and s.pathogenic_estimate=1 and s.reliability_estimate=1 "
		ret = search_db(cursor,qry)
		if not ret: continue

		patients_with_path_mutations_in_gene_2 = {}
		for p,m,s in ret:
			if not patients_with_path_mutations_in_gene_2.has_key(p): patients_with_path_mutations_in_gene_2[p] = []
			patients_with_path_mutations_in_gene_2[p].append([m,s])

		patients = 0
		p53path = 0
		miss    = 0
		cooc    = 0
		for p, muts_specs in patients_with_path_mutations_in_gene_2.iteritems():
			p53_pathogenic = False
			gene_2_missense = False
			for mutation,specimen in muts_specs:
				# is p53 pathogenic or not?
				[gist,details] = find_53_status(cursor,tumor_short,specimen)
				if gist=='pathogenic': p53_pathogenic = True
				# is mutation missense?
				qry  = "select consequence from mutations_chrom_%s " % chromosome
				qry += "where icgc_mutation_id='%s' " % mutation
				ret = search_db(cursor,qry)
				if ret:
					consequence = ret[0][0]
					if 'missense' in consequence: gene_2_missense = True

			patients += 1
			if p53_pathogenic:  p53path += 1
			if gene_2_missense: miss += 1
			if p53_pathogenic and gene_2_missense: cooc += 1

		expected = float(p53path)/patients*miss
		p_smaller, p_bigger = myfisher(patients, p53path, miss, cooc)

		total_donors  += patients
		total_p53path += p53path
		total_miss += miss
		total_cooc += cooc
		print "patients ", patients
		print "p53path  ", p53path
		print "miss     ", miss
		print "co-ocurrence:", cooc
		print "    expected: %.1f" % expected
		print "   p_smaller: %.2f" % p_smaller
		print "    p_bigger: %.2f" % p_bigger
		print
		#
		# if write_to_file: outf.write("%s\t%d\t%d\t%d\t%d\t%.1f\t%.1e\t%.1e\n"%
		#                 (tumor_short,donors, patients_with_muts_in_gene.get(gene_1, 0),
		# 				gene_2_missense, cooc,expected,p_smaller,p_bigger))
		# total_cooc += cooc

	if True:
		p_smaller, p_bigger = myfisher(total_donors, total_p53path, total_miss, total_cooc)
		print
		print "================================="
		print gene_2
		print "total donors:", total_donors
		print "tota_p53path:", total_p53path
		print "total_miss:  ", total_miss
		print "        cooc:", total_cooc
		print "    expected: %.1f" % (float(total_p53path)/total_donors*total_miss)
		print "   p_smaller: %.1e" % p_smaller
		print "    p_bigger: %.1e" % p_bigger
		expected =(float(total_p53path)/total_donors*total_miss)
		if write_to_file: outf.write("%s\t%d\t%d\t%d\t%d\t%.1f\t%.1e\t%.1e\n"%
									("total", total_donors, total_p53path,
									total_miss, total_cooc, expected, p_smaller, p_bigger))
		if write_to_file: outf.close()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
