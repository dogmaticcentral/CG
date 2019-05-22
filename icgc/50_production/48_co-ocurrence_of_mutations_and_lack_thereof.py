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
from icgc_utils.processes import  *
from config import Config
from random import shuffle

verbose = False

# ./icgc_utils/kernprof.py -l 48_....py
# python3 -m line_profiler 48_....py.lprof
# @profile
######################
def cooccurrence(tables, other_args):
	
	[bg_gene, outdir] = other_args

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#########################
	switch_to_db(cursor,"icgc")

	for table in tables:
		t0 = time.time()
		tumor_short = table.split("_")[0]

		patients_with_muts_in_gene = patients_per_gene_breakdown(cursor, table)
		total_number_of_donors = len(get_donors(cursor, table))
		nr_of_donors_w_bg_gene_mutated  = patients_with_muts_in_gene.get(bg_gene,0)
		print("="*20, "\n", table, "donors: ", total_number_of_donors, "donors with mutated %s:"%bg_gene, nr_of_donors_w_bg_gene_mutated)
		if nr_of_donors_w_bg_gene_mutated==0: continue
		# using a view here doesn't do jack squat
		# (lesson: don't trust forums)
		# bg_view = create_gene_view(cursor, table, bg_gene)
		bg_temp = create_gene_temp(cursor, table, bg_gene)
		outlines = []
		ctr = 0
		for gene in patients_with_muts_in_gene.keys():
			ctr += 1
			nr_of_donors_w_gene_mutated = patients_with_muts_in_gene.get(gene,0)
			if nr_of_donors_w_gene_mutated==0: continue
			if gene==bg_gene: continue
			if ctr%5000==0:
				print("\t\t {}: {} out of {}, {} mins".
					format(tumor_short, ctr, len(patients_with_muts_in_gene), "%.1f"%(float(time.time()-t0)/60)))
			#cooc = co_occurrence_count(cursor,table, bg_gene, gene)
			cooc = co_occurrence_count_with_a_temp(cursor, table, bg_temp, gene)
			outlines.append("%s\t%d\t%d\t%d\t%d" % (gene, total_number_of_donors,
													nr_of_donors_w_bg_gene_mutated, nr_of_donors_w_gene_mutated, cooc))
		#drop_view(cursor, bg_view)
		check_and_drop(cursor, 'icgc', bg_temp)

		if len(outlines)>0:
			outf = open("{}/{}.tsv".format(outdir,tumor_short),"w")
			outf.write("\n".join(outlines))
			outf.write("\n")
			outf.close()

		print(table, "donors: ", total_number_of_donors, "done, %.1f mins" % (float(time.time()-t0)/60))

	cursor.close()
	db.close()

#########################################
#########################################
def main():

	if len(sys.argv) < 2:
		print("usage: %s <bg gene> " % sys.argv[0])
		exit()

	bg_gene = sys.argv[1].upper()
	outdir = "cooccurrence"
	if not os.path.exists(outdir): os.mkdir(outdir)

	# divide by cancer types, because I have duplicates within each cancer type
	# that I'll resolve as I go, but I do not want the threads competing)
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in search_db(cursor,qry)]


	# ###################################
	# #number_of_chunks = 1
	#tables = ['UCEC_simple_somatic', 'BLCA_simple_somatic', 'THCA_simple_somatic', 'BRCA_simple_somatic']
	number_of_chunks = 12
	processes = parallelize(number_of_chunks, cooccurrence, tables, [bg_gene, outdir])
	if processes: wait_join(processes)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
