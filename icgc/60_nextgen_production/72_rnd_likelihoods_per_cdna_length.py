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

###############
def avg_for_random_gene_sel(cursor, table, all_genes, nr_sim_steps):

	base_qry = "select count(distinct icgc_sample_id) from %s " % table
	base_qry += "where pathogenicity_estimate=1 and reliability_estimate=1 "

	bin_population = {}
	avg_bin    = {}
	avg_sq_bin = {}

	for sel_size in range(5,310,10):
		for rep in range(nr_sim_steps):
			random_genes = sample(all_genes, sel_size)
			cdna_length  = gene_group_cdna_length(cursor, random_genes)
			bin_idx = int(cdna_length/1000)
			if not bin_idx in bin_population:
				bin_population[bin_idx] = 0
				avg_bin[bin_idx] = 0
				avg_sq_bin[bin_idx] = 0
			bin_population[bin_idx] += 1
			random_string = ",".join([quotify(g) for g in random_genes])
			qry = base_qry + "and gene_symbol in (%s)" % random_string
			mutated = error_intolerant_search(cursor, qry)[0][0]
			avg_bin[bin_idx] += mutated
			avg_sq_bin[bin_idx] += mutated**2

	return (bin_population, avg_bin, avg_sq_bin)


###############
def store_stats(cursor, tumor_short, stats_id, stats):
	bin_population, avg_bin, avg_sq_bin = stats
	for bin_idx in sorted(bin_population.keys()):
		if bin_population[bin_idx]<5: continue
		pop = bin_population[bin_idx]
		avg = avg_bin[bin_idx]/pop
		avg_sq = avg_sq_bin[bin_idx]/pop
		stdev  = sqrt(avg_sq-avg*avg)

		parameters = "{};{};{}".format(tumor_short, bin_idx*1000, pop)
		stats_string = "%.1f;%.1f"%(avg, stdev)
		#print("storing", parameters, stats_string)
		store_without_checking(cursor, 'stats',{'stats_id':stats_id, 'parameters':parameters, 'stats':stats_string})
	return

###############
def avg_pll_chunk(tables, other_args):
	[all_genes, nr_sim_steps, stats_id] = other_args
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, 'icgc')

	print("{} working on".format(os.getpid()), tables)

	for table in tables:
		tumor_short = table.split("_")[0]
		time0 = time.time()
		stats  = avg_for_random_gene_sel(cursor, table, all_genes, nr_sim_steps)
		store_stats(cursor, tumor_short, stats_id, stats)
		print("%s done in %.1f mins" %(table,  (float(time.time()-time0))/60) )

	cursor.close()
	db.close()

	return


######################################
def store_stats_description(cursor, stats_id):
	fixed_fields  = {'stats_id':stats_id}
	descr = "Random selection sample coverage per cDNA length: in how many samples will a gene from random selection of genes be mutated, " \
			"given the total cDNA length of the sample? To be compared with non-random, pathway-related selection of genes. " \
	        "Bin width: 1000 (nucleotides)."
	update_fields = {'description':descr,
					'parameters':"tumor_short:string;bin_start:int;pop_size:int",
					'stats':"average:float;stdev:float"}
	store_or_update(cursor, 'stats_description', fixed_fields, update_fields)

####################################################
def main():

	nr_sim_steps = 200

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor, 'icgc')

	# all protein coding genes
	qry = "select distinct(gene_symbol) from mutation2gene "
	all_genes = [r[0] for r in hard_landing_search(cursor, qry)]

	# how many donors/specimens have a mutation in this module
	tables = get_somatic_variant_tables(cursor)
	table_sizes = get_table_size(cursor,'icgc',tables, as_list=True)

	# we would like to store this run to our database, not to leave it laying around
	stats_id = "RSSCcdna"
	store_stats_description(cursor, stats_id)

	# random sampling
	number_of_chunks=10

	other_args = [all_genes, nr_sim_steps, stats_id]
	parallelize(number_of_chunks, avg_pll_chunk, tables, other_args, strategy='weighted', weights=table_sizes)

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
