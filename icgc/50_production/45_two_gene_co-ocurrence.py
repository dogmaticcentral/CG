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

from icgc_utils.common_queries import *
from icgc_utils.icgc_stats  import *
from config import Config
from numpy  import cumsum
from time   import time
from math import log10
verbose = False



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
		# mut_count = mutation_count_per_donor(cursor, table)
		mut_count = genes_per_patient_breakdown(cursor, table)

		patients_with_muts_in_gene = patients_per_gene_breakdown(cursor, table)
		if patients_with_muts_in_gene.get(bg_gene,0)==0: continue
		no_mutant = True
		for gene in other_genes:
			if patients_with_muts_in_gene.get(gene,0)==0: continue
			no_mutant = False
			break
		if no_mutant: continue

		cumulative_size = [0]


		total_patients = len(mut_count)
		cumulative_size.extend(cumsum(list(mut_count.values())))
		pancan_mut_count_values.extend(list(mut_count.values()))

		# hypermutated samples completely skew the stats
		# shaving off the peaks this way seems to work,
		# in the sense that is does not heavily overestimate the expected overlap
		# but it returns the results that look just like Fisher
		# should I just mark the hypermutated cancers as unreliable?
		# Jaime Iranzo, Iñigo Martincorena, and Eugene V. Koonin
		# https://www.pnas.org/content/115/26/E6010
		# consider  a sample with >3,000 mutations in a coding region a hypermutator
		# I think I should go even lower - to >1,000 mutations
		# weights = [int(10*log10(m)) for m in mut_count.values()]
		# cumulative_size.extend(cumsum(weights))
		# pancan_mut_count_values.extend(weights)

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

		expected = float(bg_gene_mutated)/total_patients*other_mutated
		print("co-ocurrence:", cooc)
		print("    expected: %.1f" % expected)
		print("   p_smaller: %.2f" % p_smaller)
		print("    p_bigger: %.2f" % p_bigger)
		if size_corrected:
			selection_sizes = [bg_gene_mutated, other_mutated, cooc]
			print("----------------------")
			p_smaller_sc, p_bigger_sc, expected_ovlp_sc = size_corrected_pvals_python(cumulative_size, selection_sizes)
			print("sc  expected (python): %.1f" % expected_ovlp_sc)
			print("sc p_smaller (python): %.2f" % p_smaller_sc)
			print("sc  p_bigger (python): %.2f" % p_bigger_sc)
			p_smaller_sc, p_bigger_sc, expected_ovlp_sc = size_corrected_pvals_C (rbf, cumulative_size, selection_sizes)
			print("sc  expected: %.1f" % expected_ovlp_sc)
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
		selection_sizes = [pancan_bg_gene, pancan_other, pancan_cooc]
		print("----------------------")
		# time0 = time()
		# python version takes 0 as the first value in the cumulative soze array
		# pancan_cumulative_size = [0]
		# pancan_cumulative_size.extend(cumsum(pancan_mut_count_values))
		# p_smaller_sc, p_bigger_sc, expected_ovlp_sc = \
		# 	size_corrected_pvals_python(pancan_cumulative_size, selection_sizes, number_of_simulation_rounds=100)
		# print("\t\t time for size corrected sim (python): %.1f mins"% (float(time()-time0)/60))
		# print("sc  expected: %.1f" % expected_ovlp_sc)
		# print("sc p_smaller: %.2e" % p_smaller_sc)
		# print("sc  p_bigger: %.2e" % p_bigger_sc)

		time0 = time()
		pancan_cumulative_size = cumsum(pancan_mut_count_values)
		p_smaller_sc, p_bigger_sc, expected_ovlp_sc = \
			size_corrected_pvals_C(rbf, pancan_cumulative_size, selection_sizes, number_of_simulation_rounds=1.0e4)
		print("\t\t time for size corrected sim: %.1f mins"% (float(time()-time0)/60))
		print("sc  expected: %.1f" % expected_ovlp_sc)
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
