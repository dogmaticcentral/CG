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


def avg_for_random_gene_sel(cursor, table, all_genes, sel_size, nr_sim_steps):
	base_qry = "select count(distinct icgc_sample_id) from %s " % table
	base_qry += "where pathogenicity_estimate=1 and reliability_estimate=1 "
	avg_random_mutated = 0
	avg_sq = 0
	for rep in range(nr_sim_steps):
		random_genes = sample(all_genes, sel_size)
		random_string = ",".join([quotify(g) for g in random_genes])
		qry = base_qry + "and gene_symbol in (%s)" % random_string
		mutated = error_intolerant_search(cursor, qry)[0][0]
		avg_random_mutated  += mutated
		avg_sq += mutated**2

	avg_random_mutated = float(avg_random_mutated)/nr_sim_steps
	avg_sq = float(avg_sq)/nr_sim_steps
	stdev = sqrt(avg_sq-avg_random_mutated*avg_random_mutated)
	return avg_random_mutated, stdev

###############
def avg_pll_chunk(tables, other_args, return_dict):
	[all_genes, selection_size, nr_sim_steps] = other_args
	avg_estimates = {}
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, 'icgc')
	for table in tables:
		avg, stdev = avg_for_random_gene_sel(cursor, table, all_genes, selection_size, nr_sim_steps)
		avg_estimates[table] = (avg,stdev)
	cursor.close()
	db.close()
	return_dict[get_process_id()]=avg_estimates


####################################################
def main():

	nr_sim_steps = 100

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor, 'icgc')

	# all protein coding genes
	qry = "select distinct(gene_symbol) from mutation2gene "
	all_genes = [r[0] for r in hard_landing_search(cursor, qry)]

	# how many donors/specimens have a mutation in this module
	tables = get_somatic_variant_tables(cursor)
	table_sizes = get_table_size(cursor,'icgc',tables, as_list=True)

	# random sampling
	outf = open("random_selection_likelihood.tsv","w")
	for selection_size in range (5,155,5):
		time0 = time.time()
		print("avg values over random gene selections of size", selection_size)
		number_of_chunks = 10
		other_args = [all_genes, selection_size, nr_sim_steps]
		avg_estimates = pll_w_return(number_of_chunks, avg_pll_chunk, tables, other_args, table_sizes)
		print("                       ... done in %.1f mins" %( (float(time.time()-time0))/60) )
		for table, stats in avg_estimates.items():
			tumor_short = table.split("_")[0]
			avg, stdev = stats
			outf.write("%d\t%s\t%.1f\t%.1f\n"%(selection_size, tumor_short,avg, stdev))
		outf.flush()

	outf.close()

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
