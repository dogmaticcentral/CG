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
import sys
from icgc_utils.processes import *
from icgc_utils.icgc_stats import *
from icgc_utils.common_queries import *
from config import Config
from math import log10
import time

def hashinit(list_of_hashes, key):
	for hash in list_of_hashes:
		hash[key] = 0

def progress_rept(nr_genes, ctr, t0):
	if (ctr%100)>0: return
	print("\t\t {}: {} out of {}, {} mins".
		format(os.getpid(), ctr, nr_genes,"%.1f" % (float(time.time() - t0) / 60)))
	return


#########################################
def gene_stats(genes, other_args):

	[rbf, precision, outdir, bg_gene, pancan_mutations, pancan_bg, pancan_cooc, pancan_mut_count_values, pancan_donors] = other_args

	nrounds=int('1'+'0'*precision)
	outf = open("{}/{}_{}_cooccurrence_probs.tsv".format(outdir, bg_gene, os.getpid()), "w")
	ctr = 0
	t0 = time.time()
	nr_genes = len(genes)
	for gene in genes:
		ctr += 1
		if pancan_mutations[gene]<10: continue
		progress_rept( nr_genes, ctr, t0)
		selection_sizes = [pancan_bg[gene], pancan_mutations[gene], pancan_cooc[gene]]
		pancan_cumulative_size = cumsum(pancan_mut_count_values[gene])
		ret = size_corrected_pvals_C(rbf, pancan_cumulative_size, selection_sizes, number_of_simulation_rounds=nrounds)
		if not ret: return None
		p_smaller, p_greater, expected_ovlp = ret

		p_smaller_fisher, p_greater_fisher = myfisher(pancan_donors[gene], pancan_bg[gene], pancan_mutations[gene], pancan_cooc[gene])
		naive_overlap = float(pancan_bg[gene])/pancan_donors[gene]*pancan_mutations[gene]
		outf.write("%s\t%d\t%d\t%d\t%d\t%.1f\t%.1e\t%.1e\t%.1f\t%.1e\t%.1e" % \
				(gene, pancan_donors[gene], pancan_mutations[gene], pancan_bg[gene],
				pancan_cooc[gene], expected_ovlp,  p_smaller, p_greater,
				 naive_overlap, p_smaller_fisher, p_greater_fisher))
		outf.write("\n")
		outf.flush()
	outf.close()

#########################################
#########################################
def main():

	if len(sys.argv) < 2:
		print("usage: %s <bg_gene>  [<precision>]" % sys.argv[0])
		print("precision - given as decimal place (3,4,5); default 2")
		exit()


	bg_gene = sys.argv[1]
	indir = "{}_cooccurrence".format(bg_gene)
	if not os.path.exists(indir) or not os.path.isdir(indir):
		print(indir, "directory not found")
		exit()

	precision = 2
	if len(sys.argv)>2:
		precision = min(int(sys.argv[2]), 6)
	# rbf is a small C program that runs the simulation
	# to evaluate Fisher-like probabilities for bins of uneven size (i.e probability of being chosen)
	rbf = Config().rbf_path()

	print("mutations per tumor per sample  ...")
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, "icgc")
	weights = {}
	qry = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in search_db(cursor, qry)]
	for table in tables:
		#mut_count[table] = mutation_count_per_donor(cursor, table)
		mut_count = genes_per_patient_breakdown(cursor, table)
		weights[table] = [int(10*log10(m)) for m in mut_count.values()]
	cursor.close()
	db.close()
	print("                                ... done")

	tsv_files = []
	for path, dir, files in  os.walk(indir):
		tsv_files += [f for f in files if f[-4:]==".tsv"]

	print("reading %s ..."%indir)
	pancan_donors = {}
	pancan_mutations = {}
	pancan_bg = {} # total bg gene mutated in tumors that have the gene of interest mutated
	pancan_cooc = {}
	pancan_mut_count_values = {}
	for file in tsv_files:
		inf = open("{}/{}".format(indir, file),"r")
		tumor_short = file.split(".")[0]
		table = "{}_simple_somatic".format(tumor_short)
		for line in inf:
			gene, donors, tot_bg, tot, tot_cooc = line.split()
			if gene not in pancan_donors:
				hashinit([pancan_donors,pancan_mutations,pancan_bg,pancan_cooc], gene)
				pancan_mut_count_values[gene] = []
			pancan_donors[gene] += int(donors)
			pancan_mutations[gene] += int(tot)
			pancan_bg[gene]        += int(tot_bg)
			pancan_cooc[gene]      += int(tot_cooc)
			#pancan_mut_count_values[gene].extend(list(mut_count[table].values()))
			pancan_mut_count_values[gene].extend(list(weights[table]))
		inf.close()
	print("            ... done")

	outdir = "{}_coocc_stats".format(bg_gene)
	if not os.path.exists(outdir): os.mkdir(outdir)
	other_args = [rbf, precision, outdir, bg_gene, pancan_mutations, pancan_bg, pancan_cooc, pancan_mut_count_values, pancan_donors]
	number_of_chunks = 10
	parallelize(number_of_chunks, gene_stats, list(pancan_mutations.keys()), other_args)



#########################################
if __name__ == '__main__':
	main()


#########################################
'''
further sorting (with rank number) 7: anticorrelate, 8: correlates

cat *.tsv | sort -gk8  | awk '{ct +=1; printf "%6d  ", ct; print}' |  grep RPL | grep -v MRPL
'''
