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


import os
from scipy import stats


def hashinit(list_of_hashes, key):
	for hash in list_of_hashes:
		hash[key] = 0


#########################################
def myfisher(donors, tp53_mutated, other_mutated, cooc):
	# https://en.wikipedia.org/wiki/Hypergeometric_distribution
	# https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html
	# In probability theory and statistics, the cumulative distribution function
	# (CDF, also cumulative density function) of a real-valued random variable X, or just distribution function of X,
	#  evaluated at x,
	# is the probability that X will take a value less than or equal to x
	p_smaller = stats.hypergeom.cdf(cooc, donors, tp53_mutated, other_mutated)
	if cooc>0: # the less or equal in both cases is the reason why they do not add up to 1
		p_greater = 1 -  stats.hypergeom.cdf(cooc-1, donors, tp53_mutated, other_mutated)
	else:
		p_greater = 1.0
	return p_smaller, p_greater

#########################################
#########################################
def main():

	if not os.path.exists("coocurrence"):
		print "coocurrence dir not found";
		exit()

	tsv_files = []
	for path, dir, files in  os.walk("coocurrence"):
		tsv_files += [f for f in files if f[-4:]==".tsv"]

	total_donors = {}
	total = {}
	total_tp53 = {} # total tp53 mutated in tumors that have gene of interest mutated
	total_cooc = {}
	for file in tsv_files:
		inf = open("coocurrence/"+file,"r")
		for line in inf:
			gene, donors, tot_tp53, tot, tot_cooc = line.split()
			if not total_donors.has_key(gene): hashinit([total_donors,total,total_tp53,total_cooc], gene)
			total_donors[gene] += int(donors)
			total[gene]        += int(tot)
			total_tp53[gene]   += int(tot_tp53)
			total_cooc[gene]   += int(tot_cooc)
		inf.close()

	ct = 0
	for gene in total.keys():
		fract_tp53 = float(total_tp53[gene])/total_donors[gene]
		expected   = total[gene]*fract_tp53
		if total[gene]<10: continue
		if expected==0: continue
		ratio = total_cooc[gene]/(expected)
		p_smaller, p_greater = myfisher(total_donors[gene], total_tp53[gene], total[gene], total_cooc[gene])
		#print "%10s   %5d %5d %5d    %5d    %6.2f   %.2f     %.1e    %.1e " % \
		#		(gene, total_donors[gene], total[gene], total_tp53[gene],
		#		total_cooc[gene], expected, ratio, p_smaller, p_greater)
		print "%s\t%d\t%d\t%d\t%d\t%.1f\t%.1e\t%.1e" % \
				(gene, total_donors[gene], total[gene], total_tp53[gene],
				total_cooc[gene], expected,  p_smaller, p_greater)
		ct += 1
		#if ct==100: break
#########################################
'''
further sorting (with rank number)
sort -gk7 anticorrelates.txt | awk '{ct +=1; printf "%6d", ct; print}' |  grep RPL | grep -v MRPL
'''

#########################################
if __name__ == '__main__':
	main()
