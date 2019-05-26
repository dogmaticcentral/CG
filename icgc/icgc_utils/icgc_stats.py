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
import os
import subprocess

from scipy  import stats
from random import randint
from numpy  import searchsorted, cumsum
from math import exp

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
def size_corrected_pvals_C(rbf, cumulative_bin_size, selection_sizes, number_of_simulation_rounds=1000):

	[selection_size_1, selection_size_2, overlap_size] = selection_sizes
	outf_name = "bdries_{}.txt".format(os.getpid())
	boutf = open(outf_name,"w")
	boutf.write("\n".join([str(i) for i in cumulative_bin_size]) + "\n")
	boutf.close()

	cmd = "{} {} {} {} {} {}".format(rbf, outf_name, selection_size_1, selection_size_2, overlap_size, number_of_simulation_rounds)
	try:
		line = subprocess.check_output(cmd, shell=True).decode('utf8').split("\n")[0]
		token = line.rstrip().split("\t")
		if token[0] != "OK":
			print("error running:\n"+cmd)
			print("returned:\n"+line)
			return None
		[pval_lt, pval_gt, expected_overlap] = [float(t) for t in token[1:4]]
		os.remove(outf_name)
		return pval_lt, pval_gt, expected_overlap

	except subprocess.CalledProcessError as xcptn:
		print("error running:\n"+cmd)
		print("returned:\n{}".format(xcptn.output))
		return None


#########################
def size_corrected_pvals_python(cumulative_size, selection_sizes, number_of_simulation_rounds=100):
	# python version takes 0 as the first value in the cumulative size array
	# for example
	# cumulative_size = [0]
	# cumulative_size.extend(cumsum(mut_count_values))
	[selection_size_1, selection_size_2, overlap_size] = selection_sizes
	count_smaller_or_equal = 0;
	count_bigger_or_equal = 0
	exptected_ovlp = 0
	for r in  range(number_of_simulation_rounds): # simulation replicates
		selection_1 = bin_selection(cumulative_size, selection_size_1)
		selection_2 = bin_selection(cumulative_size, selection_size_2)
		random_ovlp_size = len(selection_1.intersection(selection_2))
		if random_ovlp_size<=overlap_size: count_smaller_or_equal+=1
		if random_ovlp_size>=overlap_size: count_bigger_or_equal+=1
		exptected_ovlp += random_ovlp_size

	pval_lt  = float(count_smaller_or_equal)/number_of_simulation_rounds
	pval_gt  = float(count_bigger_or_equal)/number_of_simulation_rounds
	ovlp     = float(exptected_ovlp)/number_of_simulation_rounds

	return pval_lt, pval_gt, ovlp



###################################
def main():
	#sample_size , selection_size_1, selection_size_2, overlap_size = [100, 50, 10, 3]
	sample_size , selection_size_1, selection_size_2, overlap_size = [1000, 500, 100, 50]
	selection_sizes = [selection_size_1, selection_size_2, overlap_size]
	p_smaller_fisher, p_bigger_fisher = fisher(sample_size, selection_size_1, selection_size_2, overlap_size)

	sample_weights = [10]*sample_size
	weights_cumulative = cumsum(sample_weights)
	#print(weights_cumulative)
	p_smaller_py, p_bigger_py, ovlp_py = size_corrected_pvals_python(weights_cumulative, selection_sizes, number_of_simulation_rounds=1000)
	p_smaller_C, p_bigger_C, ovlp_C = size_corrected_pvals_python(weights_cumulative, selection_sizes, number_of_simulation_rounds=1000)

	print("Fisher: %.2e  %.2e"%(p_smaller_fisher, p_bigger_fisher))
	print("Python: %.2e  %.2e"%(p_smaller_py, p_bigger_py))
	print("     C: %.2e  %.2e"%(p_smaller_C, p_bigger_C))


	sample_weights = [int(sample_weights[i]*exp(-abs(i-sample_size/2)/500))for i in range(sample_size)]
	weights_cumulative = cumsum(sample_weights)

	p_smaller_py, p_bigger_py, ovlp_py = size_corrected_pvals_python(weights_cumulative, selection_sizes, number_of_simulation_rounds=1000)
	p_smaller_C, p_bigger_C, ovlp_C = size_corrected_pvals_python(weights_cumulative, selection_sizes, number_of_simulation_rounds=1000)

	print("Fisher: %.2e  %.2e"%(p_smaller_fisher, p_bigger_fisher))
	print("Python: %.2e  %.2e"%(p_smaller_py, p_bigger_py))
	print("     C: %.2e  %.2e"%(p_smaller_C, p_bigger_C))






#########################################
if __name__ == '__main__':
	main()

