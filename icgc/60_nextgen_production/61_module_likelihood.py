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
import time

from icgc_utils.common_queries import *
from icgc_utils.processes import *
from config import Config
from random import sample
from math import sqrt

def gulp_in_list(listname, column=0):
	if not os.path.exists(listname):
		print(listname, "not found")
		exit()
	with open(listname,"r") as inf:
		retlist = [line.split("\t")[column].strip() for line in inf]
	return retlist

def avg_for_random_gene_sel(cursor, table, all_genes, sel_size):
	base_qry = "select count(distinct icgc_sample_id) from %s " % table
	base_qry += "where pathogenicity_estimate=1 and reliability_estimate=1 "
	avg_random_mutated = 0
	avg_sq = 0
	for rep in range(100):
		random_genes = sample(all_genes, sel_size)
		random_string = ",".join([quotify(g) for g in random_genes])
		qry = base_qry + "and gene_symbol in (%s)" % random_string
		mutated = error_intolerant_search(cursor, qry)[0][0]
		avg_random_mutated  += mutated
		avg_sq += mutated**2

	avg_random_mutated = float(avg_random_mutated)/100
	avg_sq = float(avg_sq)/100
	stdev = sqrt(avg_sq-avg_random_mutated*avg_random_mutated)
	return avg_random_mutated, stdev

###############
def avg_pll_chunk(tables, other_args, return_dict):
	[all_genes, selection_size] = other_args
	avg_estimates = {}
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, 'icgc')
	for table in tables:
		avg, stdev = avg_for_random_gene_sel(cursor, table, all_genes, selection_size)
		avg_estimates[table] = (avg,stdev)
	cursor.close()
	db.close()
	return_dict[get_process_id()]=avg_estimates


####################################################
def main():

	if len(sys.argv) < 2:
		print("usage: %s <input list, tsv> [<column>]" % sys.argv[0])
		exit()
	infilename = sys.argv[1]
	gene_list = gulp_in_list(infilename, int(sys.argv[2])-1 if len(sys.argv)>2 else 0)
	print("number of genes in {}: {}".format(infilename, len(gene_list)))
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	# check input sanity
	for gene in gene_list:
		if check_gene_hgnc(cursor, gene): continue
		print(gene, "not found among HGNC approved symbols table")
		exit()

	switch_to_db(cursor, 'icgc')

	# which of our input genes are ever mutated?
	mutated_gene_list = set(gene_list)

	for gene in gene_list:
		qry = "select * from mutation2gene where gene_symbol='%s' limit 1" % gene
		ret = error_intolerant_search(cursor, qry)
		if ret: continue
		print(gene, "never mutated")
		mutated_gene_list.remove(gene)

	# all protein coding genes
	qry = "select distinct(gene_symbol) from mutation2gene "
	all_genes = [r[0] for r in hard_landing_search(cursor, qry)]

	# how many donors/specimens have a mutation in this module
	tables = get_somatic_variant_tables(cursor)
	table_sizes = get_table_size(cursor,'icgc',tables, as_list=True)

	# random sampling, for comparison
	time0 = time.time()
	print("avg values over random gene selections of the same size ... ")
	number_of_chunks = 10
	other_args = [all_genes, len(mutated_gene_list)]
	avg_estimates = pll_w_return(number_of_chunks, avg_pll_chunk, tables, other_args, table_sizes)
	print("                       ... done in %.1f mins" %( (float(time.time()-time0))/60) )

	gene_string = ",".join([quotify(g) for g in mutated_gene_list])
	for table in tables:
		base_qry = "select count(distinct icgc_sample_id) from %s " % table
		base_qry += "where pathogenicity_estimate=1 and reliability_estimate=1 "
		total_donor_in_table = error_intolerant_search(cursor, base_qry)[0][0]

		qry = base_qry + "and gene_symbol in (%s)" % gene_string
		module_mutated  = error_intolerant_search(cursor, qry)[0][0]

		avg_random_mutated, stdev = avg_estimates[table]
		print("{}: {} donors".format(table,  total_donor_in_table))
		z = (module_mutated-avg_random_mutated)/stdev
		print("            %5d   %5.1f  %5.2f"% (module_mutated, avg_random_mutated, z))



	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
