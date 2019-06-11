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

from config import Config
from icgc_utils.reactome import *
from icgc_utils.common_queries import *
from icgc_utils.utils import cancer_dictionary
import re
####################################################
def parse(infile):
	gene_groups = []
	inf = open(infile, "r")
	entry = []
	for line in inf:
		line = line.strip()
		if 'expected' in line:
			name = entry[0]
			[number_of_genes, cdna_length] = map(int, re.findall(r'\d+', entry[1]))
			[donors_affected] = map(int, re.findall(r'\d+', entry[2]))
			[expected_donors, stdev, zscore] =  map(float, re.findall(r'\d+\.\d+', line))
			gene_groups.append([name, number_of_genes, cdna_length, donors_affected, expected_donors, stdev, zscore])
		elif len(line)==0:
			entry = []
		else:
			entry.append(line)

	gene_groups.sort(key=lambda x:x[-1], reverse=True)
	return gene_groups

####################################################
def elaborate (cursor, table, reactome_gene_groups, name2reactome_id, group):
	[name, number_of_genes, cdna_length, donors_affected, expected_donors, stdev, zscore] = group
	if zscore<5: return
	reactome_id = name2reactome_id[name]
	print(name, zscore)
	print("\t number of genes:  %3d   combined cdna length: %d " % (number_of_genes, cdna_length) )
	print("\t number of donors affected:  %3d " %  donors_affected)
	print("\t expected donors affected:  %6.2f      stdev:  %6.2f    z:  %6.2f " % (expected_donors, stdev, zscore))
	if donors_affected==0: return

	print("all genes", len(reactome_gene_groups[reactome_id]))
	print(sorted(reactome_gene_groups[reactome_id]))
	gene_string = ",".join([quotify(g) for g in reactome_gene_groups[reactome_id]])
	qry  = "select distinct(gene_symbol) from %s " % table
	qry += "where pathogenicity_estimate=1 and reliability_estimate=1 "
	qry += "and gene_symbol in (%s)" % gene_string
	genes_mutated = [r[0] for r in error_intolerant_search(cursor, qry)]
	print("actually affected", len(genes_mutated))
	print(sorted(genes_mutated))

	print()

####################################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, 'icgc')

	cancer_dict = cancer_dictionary()
	reactome_gene_groups = find_gene_groups(cursor)
	name2reactome_id = dict(hard_landing_search(cursor, "select name, reactome_pathway_id from reactome_pathways"))

	indir = "gene_groups"
	for tumor_short in sorted(os.listdir(indir)):
		table = tumor_short+"_simple_somatic"
		print("\n=============================")
		print(tumor_short)
		if tumor_short in cancer_dict: print(cancer_dict[tumor_short]["description"])
		print("total donors:", len(get_donors(cursor, table)))
		print()
		gene_groups = parse("{}/{}/{}".format(indir, tumor_short, 'pathways.txt'))
		for group in gene_groups:
			elaborate(cursor, table, reactome_gene_groups, name2reactome_id,  group)

	cursor.close()
	db.close()

	return


#########################################
if __name__ == '__main__':
	main()
