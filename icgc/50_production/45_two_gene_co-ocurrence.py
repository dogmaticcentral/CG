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
from scipy  import stats
from config import Config
from random import randint
from numpy  import searchsorted, cumsum
from time   import time
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

#####################################
def bin_selection(cumulative_size, number_of_selections):
	selection = set([])
	for b in range(number_of_selections):
		while True: # no replacement - otherwise the most probable bins are the attractors
			# random.randint(a, b) Return a random integer N such that a <= N <= b
			random_val = randint(1, cumulative_size[-1])
			# numpy.searchsorted(a, v, side='left', sorter=None)
			# side='left' : left 	a[i-1] < v <= a[i]
			bin_this_val_belongs_to = searchsorted(cumulative_size, random_val)
			if not bin_this_val_belongs_to in selection: break
		selection.add(bin_this_val_belongs_to)
	return selection

#####################################
def size_corrected_pvals_C(rbf, cumulative_size, selection_size_1, selection_size_2, overlap_size):

	number_of_simulation_rounds = 1000
	pval_lt, pval_gt = 0, 0
	boutf = open("bdries.txt","w")
	boutf.write("\n".join([str(i) for i in cumulative_size]) + "\n")
	boutf.close()

	cmd = "{} bdries.txt {} {} {} {}".format(rbf, selection_size_1, selection_size_2, overlap_size, number_of_simulation_rounds)
	line = subprocess.check_output(cmd, shell=True).decode('utf8').split("\n")[0]
	token = line.rstrip().split("\t")
	if token[0] != "OK":
		print("error running:\ncmd")
		print("returned:\nline")
		exit()
	[pval_lt, pval_gt] = [float(t) for t in token[1:3]]
	return pval_lt, pval_gt


def size_corrected_pvals_python(cumulative_size, selection_size_1, selection_size_2, overlap_size):
	#print(cumulative_size, bg_gene_mutated, other_mutated, cooc)
	number_of_simulation_rounds = 100
	count_smaller_or_equal = 0;
	count_bigger_or_equal = 0
	for r in  range(number_of_simulation_rounds): # simulation replicates
		selection_1 = bin_selection(cumulative_size, selection_size_1)
		selection_2 = bin_selection(cumulative_size, selection_size_2)
		random_ovlp_size = len(selection_1.intersection(selection_2))
		if random_ovlp_size<=overlap_size: count_smaller_or_equal+=1
		if random_ovlp_size>=overlap_size: count_bigger_or_equal+=1

	pval_lt  = float(count_smaller_or_equal)/number_of_simulation_rounds
	pval_gt  = float(count_bigger_or_equal)/number_of_simulation_rounds

	return pval_lt, pval_gt

###################################
def main():

	size_corrected = True

	if len(sys.argv) < 3:
		print("usage: %s <bg gene>  <gene 1> [<gene 2> ...]" % sys.argv[0])
		exit()

	bg_gene = sys.argv[1].upper()
	other_genes = [g.upper() for g in sys.argv[2:]]

	# TODO: check that all gene names exist

	# rbf is a small C program that runs the simulation
	# to evaluate Fisher-like probabilities for bins of uneven size (i.e probaility of being chosen)

	rbf = Config().rbf_path()


	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]

	for table in tables:
		#print("checking/creating index on", table)
		create_index(cursor, 'icgc', 'donor_gene_idx', table,['icgc_donor_id','gene_symbol'])
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

	pancan_mut_count_values = []
	p_smaller_sc, p_bigger_sc = 0 ,0 # to make the code checker shut up
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

		cumulative_size = [0]
		mut_count = mutation_count_per_donor(cursor, table)
		total_patients = len(mut_count)
		cumulative_size.extend(cumsum(list(mut_count.values())))
		pancan_mut_count_values.extend(list(mut_count.values()))

		print("=================================")
		print(table)
		print("donors: ", total_patients)

		pancan_donors += total_patients
		for gene in [bg_gene]+other_genes:
			print(gene, patients_with_muts_in_gene.get(gene, 0))
		#cumulative_size.append(cumulative_size[-1]+total_patients)

		bg_gene_mutated  = patients_with_muts_in_gene.get(bg_gene,0)
		other_mutated    = patients_with_muts_in_gene_group(cursor, table, other_genes)
		pancan_bg_gene  += bg_gene_mutated
		pancan_other    += other_mutated

		cooc = co_ocurrence_w_group_count(cursor, table, bg_gene, other_genes)

		p_smaller, p_bigger = myfisher(total_patients, bg_gene_mutated, other_mutated, cooc)
		#pval_lt, pval_gt = fisher(donors, gene_1_mutated, other_mutated, cooc)
		if size_corrected: p_smaller_sc, p_bigger_sc = size_corrected_pvals_C (rbf, cumulative_size, bg_gene_mutated, other_mutated, cooc)

		expected = float(bg_gene_mutated)/total_patients*other_mutated
		print("co-ocurrence:", cooc)
		print("    expected: %.1f" % expected)
		print("   p_smaller: %.2f" % p_smaller)
		print("    p_bigger: %.2f" % p_bigger)
		if size_corrected:
			print("----------------------")
			print("sc p_smaller: %.2f" % p_smaller_sc)
			print("sc  p_bigger: %.2f" % p_bigger_sc)

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
	if size_corrected:
		time0 = time()
		# python version takes 0 as the first value in the cumulative soze array
		# pancan_cumulative_size = [0]
		# pancan_cumulative_size.extend(cumsum(pancan_mut_count_values))
		# p_smaller_sc, p_bigger_sc = size_corrected_pvals_python(pancan_cumulative_size, pancan_bg_gene, pancan_other, pancan_cooc)
		pancan_cumulative_size = cumsum(pancan_mut_count_values)
		p_smaller_sc, p_bigger_sc = size_corrected_pvals_C(rbf, pancan_cumulative_size, pancan_bg_gene, pancan_other, pancan_cooc)
		print("----------------------")
		print("\t\t time for szie corrected sim: %.1f mins"% (float(time()-time0)/60))
		print("sc p_smaller: %.2e" % p_smaller_sc)
		print("sc  p_bigger: %.2e" % p_bigger_sc)
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
