#! /usr/bin/python3
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
import subprocess

from icgc_utils.common_queries import *
from scipy import stats
from config import  Config

verbose = False


#####################################
def myfisher(donors, gene_1_mutated, other_mutated, cooc):
	# https://en.wikipedia.org/wiki/Hypergeometric_distribution
	# https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html
	# In probability theory and statistics, the cumulative distribution function
	# (CDF, also cumulative density function) of a real-valued random variable X, or just distribution function of X,
	#  evaluated at x,
	# is the probability that X will take a value less than or equal to x
	p_smaller = stats.hypergeom.cdf(cooc, donors, gene_1_mutated, other_mutated)
	if cooc>0: # the less or equal in both cases is the reason why they do not add up to 1
		p_greater = 1 - stats.hypergeom.cdf(cooc-1, donors, gene_1_mutated, other_mutated)
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
def main():

	if len(sys.argv) < 3:
		print("usage: %s <bg gene>  <gene 1> [<gene 2> ...]" % sys.argv[0])
		exit()

	bg_gene = sys.argv[1].upper()
	other_genes = [g.upper() for g in sys.argv[2:]]

	# TODO: check that all gene names exist

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]

	#########################
	switch_to_db(cursor,"icgc")
	pancan_donors = 0
	pancan_other  = 0
	pancan_bg_gene = 0
	pancan_cooc = 0

	write_to_file =  len(other_genes)==1

	if write_to_file:
		outf = open("{}_cooccurrence.tsv".format("_".join([bg_gene]+other_genes)),"w")
		outf.write("\t".join(['cancer','donors', "mutations in %s"%bg_gene,
							"mutations in %s"%other_genes[0], 'cooccurrence','expected',
							'p smaller', 'p bigger'])+"\n")

	for table in tables:
		tumor_short = table.split("_")[0]
		patients_with_muts_in_gene = patients_per_gene_breakdown(cursor, table)
		if patients_with_muts_in_gene.get(bg_gene,0)==0: continue
		no_mutant = True
		for gene in other_genes:
			if patients_with_muts_in_gene.get(gene,0)==0: continue
			no_mutant = False
			break
		if no_mutant: continue

		total_patients = len(get_donors(cursor, table))

		print("=================================")
		print(table)
		print("donors: ", total_patients)

		pancan_donors += total_patients
		for gene in [bg_gene]+other_genes:
			print(gene, patients_with_muts_in_gene.get(gene, 0))

		bg_gene_mutated  = patients_with_muts_in_gene.get(bg_gene,0)
		other_mutated    = patients_with_muts_in_gene_group(cursor, table, other_genes)
		pancan_bg_gene  += bg_gene_mutated
		pancan_other    += other_mutated

		cooc = co_ocurrence_w_group_count(cursor, table, bg_gene, other_genes)

		p_smaller, p_bigger = myfisher(total_patients, bg_gene_mutated, other_mutated, cooc)
		#pval_lt, pval_gt = fisher(donors, gene_1_mutated, other_mutated, cooc)

		expected = float(bg_gene_mutated)/total_patients*other_mutated
		print("co-ocurrence:", cooc)
		print("    expected: %.1f" % expected)
		print("   p_smaller: %.2f" % p_smaller)
		print("    p_bigger: %.2f" % p_bigger)
		print()

		if write_to_file: outf.write("%s\t%d\t%d\t%d\t%d\t%.1f\t%.2f\t%.1f\n"%
						(tumor_short,total_patients, patients_with_muts_in_gene.get(bg_gene, 0),
						patients_with_muts_in_gene.get(other_genes[0], 0),
						cooc,expected,p_smaller,p_bigger))
		pancan_cooc += cooc


	p_smaller, p_bigger = myfisher(pancan_donors, pancan_bg_gene, pancan_other, pancan_cooc)
	print()
	print("=================================")
	print(other_genes)
	print("total donors:", pancan_donors)
	print("       other:", pancan_other)
	print("%12s: %d" % (bg_gene, pancan_bg_gene))
	print("        cooc:", pancan_cooc)
	print("    expected: %.1f" % (float(pancan_bg_gene)/pancan_donors*pancan_other))
	print("   p_smaller: %.1e" % p_smaller)
	print("    p_bigger: %.1e" % p_bigger)
	expected = (float(pancan_bg_gene)/pancan_donors*pancan_other)
	if write_to_file: outf.write("%s\t%d\t%d\t%d\t%d\t%.1f\t%.1e\t%.1e\n"%
								("total", pancan_donors, pancan_bg_gene,
								pancan_other, pancan_cooc, expected, p_smaller, p_bigger))


	if write_to_file: outf.close()

	#print myfisher(total_donors*4, total_gene_1*4, total_other*4, total_cooc*4)
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
