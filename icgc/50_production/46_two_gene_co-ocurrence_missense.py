#! /usr/bin/python
#
# This source code is part of icgc, an ICGC processing pipeline.
# 
# Icgc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Icgc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#

# coocurrence of missense (only)  with wildtype

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
def patients_with_missense(cursor, somatic_table, gene_2):

	mutations_table = "mutations_chrom_%s" % find_chromosome(cursor, gene_2)
	qry  = "select count(distinct(s.icgc_donor_id)) from %s s, %s m, mutation2gene g " % (somatic_table, mutations_table)
	qry += "where g.gene_symbol='%s' " % gene_2
	qry += "and g.icgc_mutation_id=s.icgc_mutation_id "
	qry += "and g.icgc_mutation_id=m.icgc_mutation_id "
	qry += "and m.consequence like '%missense%' "
	qry += "and s.reliability_estimate=1 "
	ret = search_db(cursor,qry)

	if ret: return ret[0][0]
	return 0

###################################
def co_ocurrence_any_vs_missense(cursor, somatic_table, gene_1, gene_2):

	mutations_table = "mutations_chrom_%s" % find_chromosome(cursor, gene_2)
	qry  = "select count(distinct(s1.icgc_donor_id)) "
	qry += "from %s s1, %s s2, %s m, mutation2gene g1, mutation2gene g2  " % (somatic_table, somatic_table, mutations_table)
	qry += "where g1.gene_symbol='%s' " % gene_1
	qry += "and g2.gene_symbol='%s' " % gene_2
	qry += "and g2.icgc_mutation_id=s2.icgc_mutation_id "
	qry += "and g2.icgc_mutation_id=m.icgc_mutation_id "
	qry += "and m.consequence like '%missense%' "
	qry += "and g1.icgc_mutation_id=s1.icgc_mutation_id "
	qry += "and s1.pathogenic_estimate=1 "
	qry += "and s1.reliability_estimate=1 "
	qry += "and s2.reliability_estimate=1 "
	qry += "and s1.icgc_donor_id=s2.icgc_donor_id "
	ret = search_db(cursor,qry)

	if ret: return ret[0][0]
	return 0

###################################
# ## file:///home/ivana/Dropbox/Sinisa/ribosomal/html/the_curious_case_of_rpl22.html
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	gene_1 = 'TP53'
	#gene_1 = 'BRCA1'
	gene_2  = 'RPL11'
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

	write_to_file =  False

	if write_to_file:
		outf = open("{}_{}_cooccurrence.tsv".format(gene_1, gene_2),"w")
		outf.write("\t".join(['cancer','donors', "mutations in %s"%gene_1,
			                "mutations in %s"%gene_2, 'cooccurrence','expected',
			                'p smaller', 'p bigger'])+"\n")

	for table in tables:
		tumor_short = table.split("_")[0]
		patients_with_muts_in_gene = patients_per_gene_breakdown(cursor, table)

		if patients_with_muts_in_gene.get(gene_1,0)==0: continue
		if patients_with_muts_in_gene.get(gene_2,0)==0: continue


		print("=================================")
		print(table)
		donors = len(get_donors(cursor, table))
		print("donors: ", donors)
		total_donors += donors
		print(gene_2, patients_with_muts_in_gene.get(gene_2, 0))

		gene_1_mutated  = patients_with_muts_in_gene.get(gene_1,0)
		gene_2_missense = patients_with_missense(cursor, table, gene_2)
		print("patients with missense:", gene_2_missense)
		total_gene_1   += gene_1_mutated
		total_other    += gene_2_missense

		cooc = co_ocurrence_any_vs_missense(cursor, table, gene_1, gene_2)

		p_smaller, p_bigger = myfisher(donors, gene_1_mutated, gene_2_missense, cooc)

		expected = float(gene_1_mutated)/donors*gene_2_missense
		print("co-ocurrence:", cooc)
		print("    expected: %.1f" % expected)
		print("   p_smaller: %.2f" % p_smaller)
		print("    p_bigger: %.2f" % p_bigger)
		print()

		if write_to_file: outf.write("%s\t%d\t%d\t%d\t%d\t%.1f\t%.1e\t%.1e\n"%
		                (tumor_short,donors, patients_with_muts_in_gene.get(gene_1, 0),
						gene_2_missense, cooc,expected,p_smaller,p_bigger))
		total_cooc += cooc


	p_smaller, p_bigger = myfisher(total_donors, total_gene_1, total_other, total_cooc)
	print()
	print("=================================")
	print(gene_2)
	print("total donors:", total_donors)
	print("        other:", total_other)
	print("          %s: %d" % (gene_1, total_gene_1))
	print("        cooc:", total_cooc)
	print("    expected: %.1f" % (float(total_gene_1)/total_donors*total_other))
	print("   p_smaller: %.1e" % p_smaller)
	print("    p_bigger: %.1e" % p_bigger)
	expected = (float(total_gene_1)/total_donors*total_other)
	if write_to_file: outf.write("%s\t%d\t%d\t%d\t%d\t%.1f\t%.1e\t%.1e\n"%
								("total", total_donors, total_gene_1,
								total_other, total_cooc, expected, p_smaller, p_bigger))


	if write_to_file: outf.close()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
