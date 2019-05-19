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


import time
from icgc_utils.common_queries import *
from icgc_utils.processes import  *
verbose = False

def hashinit(list_of_hashes, key):
	for hash in list_of_hashes:
		hash[key] = 0


######################
def coocurrence(tables, other_args):

	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	switch_to_db(cursor,"icgc")

	for table in tables:
		t0 = time.time()
		tumor_short = table.split("_")[0]
		patients_with_muts_in_gene = patients_per_gene_breakdown(cursor, table)
		donors = len(get_donors(cursor, table))
		print("="*20, "\n", table, "donors: ", donors, "donors with mutated TP53:", patients_with_muts_in_gene.get('TP53',0))
		if patients_with_muts_in_gene.get('TP53',0)==0: continue

		# the total sums are per gene, over tumors that have a mutation in this gene
		total_donors = {}
		total = {}
		total_tp53 = {} # total tp53 mutated in tumors that have gene of interest mutated
		total_cooc = {}
		for gene, number_of_patients in patients_with_muts_in_gene.items():
			if gene=='TP53': continue
			if gene not in total_donors: hashinit([total_donors,total,total_tp53,total_cooc], gene)
			total_donors[gene] += donors
			total_tp53[gene]   += patients_with_muts_in_gene.get('TP53',0)
			total[gene] += number_of_patients
			total_cooc[gene] += co_ocurrence_count(cursor, table, gene, 'TP53')

		print(table, "donors: ", donors, "done, %.1f mins" % (float(time.time()-t0)/60))

		outf = open("coocurrence/%s.tsv"%tumor_short,"w")
		for gene, donors in total_donors.items():
			if donors<10: continue
			outf.write("%s\t%d\t%d\t%d\t%d\n" % (gene, donors, total_tp53[gene], total[gene], total_cooc[gene]))
		outf.close()

	cursor.close()
	db.close()

#########################################
#########################################
def main():

	if not os.path.exists("coocurrence"): os.mkdir("coocurrence")

	# divide by cancer types, because I have duplicates within each cancer type
	# that I'll resolve as I go, but I do not want the threads competing)
	db     = connect_to_mysql()
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in search_db(cursor,qry)]

	cursor.close()
	db.close()

	#tables= ['ALL_simple_somatic']
	number_of_chunks = 8  # myISAM does not deadlock
	parallelize(number_of_chunks, coocurrence, tables, [])



#########################################
if __name__ == '__main__':
	main()
