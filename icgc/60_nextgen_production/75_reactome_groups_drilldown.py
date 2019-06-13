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
			[expected_donors, stdev, zscore] =  map(float, re.findall(r'-*\d+\.\d+', line))
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
	if abs(zscore)<5: return None
	reactome_id = name2reactome_id[name]
	retlines = []
	retlines.append("{}  {}".format(name, zscore))
	retlines.append("\t number of genes:  %3d   combined cdna length: %d " % (number_of_genes, cdna_length) )
	retlines.append("\t number of donors affected:  %3d " %  donors_affected)
	retlines.append("\t expected donors affected:  %6.2f      stdev:  %6.2f    z:  %6.2f " % (expected_donors, stdev, zscore))
	if donors_affected==0: return

	retlines.append("all genes in the Reactome group %d " % len(reactome_gene_groups[reactome_id]))
	retlines.append(str(sorted(reactome_gene_groups[reactome_id])))
	gene_string = ",".join([quotify(g) for g in reactome_gene_groups[reactome_id]])
	qry  = "select  gene_symbol, count(distinct(icgc_donor_id)) from %s " % table
	qry += "where pathogenicity_estimate=1 and reliability_estimate=1 "
	qry += "and gene_symbol in (%s) group by gene_symbol" % gene_string
	# Under Python 3.6, the built-in dict does track insertion order,
	# although this behavior is a side-effect of an implementation change and should not be relied on.
	genes_mutated = dict(sorted(error_intolerant_search(cursor, qry), key= lambda r: r[1], reverse=True))
	retlines.append("actually affected %d " % len(genes_mutated))
	retlines.append(str(genes_mutated))
	retlines.append("")
	return retlines

###################
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def count_silent(codon):
	orig_translation = Seq(codon, generic_dna).translate()
	orig_codon = list(codon)
	silent = 0
	nonsilent = 0
	for i in range(len(orig_codon)):
		for substitution in ['A', 'C', 'T', 'G']:
			if substitution == orig_codon[i]: continue
			mutated_codon = orig_codon.copy()
			mutated_codon[i] = substitution
			translation = Seq(''.join(mutated_codon), generic_dna).translate()
			if translation==orig_translation:
				silent+=1
			else:
				nonsilent += 1
	return silent, nonsilent

####################################################
def silent_nonsilent_theoretical(cursor, gene_symbol):
	ratio = 1
	gene_symbol = gene_symbol.upper()
	qry = "select s.sequence from ensembl_coding_seqs s, ensembl_ids i, hgnc h "
	qry += "where s.transcript_id=i.transcript and i.gene=h.ensembl_gene_id and h.approved_symbol='%s' " % gene_symbol
	ret = error_intolerant_search(cursor,qry)
	if not ret: return -1
	sequence = ret[0][0]
	print(len(sequence)/3)
	read_in_last_codon = False
	silent = 0
	nonsilent = 0
	for codon in [sequence[i:i+3] for i in range(0, len(sequence),3)]:
		if read_in_last_codon:
			print("{} the sequence  has a nonsense stop codon".format(gene_symbol))
			exit()
		bp_codon = Seq(codon, generic_dna)
		translation = bp_codon.translate()
		read_in_last_codon = translation=="*"
		[codon_silent, codon_nonsilent] = count_silent(codon)
		print(codon, translation, codon_silent, codon_nonsilent)
		silent += codon_silent
		nonsilent += codon_nonsilent
	return silent, nonsilent

####################################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, 'icgc')

	silent, nonsilent = silent_nonsilent_theoretical(cursor, '')
	print(" %d   %d  %.3f   %.3f"%(silent, nonsilent, silent/nonsilent, silent/(silent+nonsilent)))
	exit()

	cancer_dict = cancer_dictionary()
	reactome_gene_groups = find_gene_groups(cursor)
	name2reactome_id = dict(hard_landing_search(cursor, "select name, reactome_pathway_id from reactome_pathways"))

	indir = "gene_groups"
	for tumor_short in sorted(os.listdir(indir)):
		print(tumor_short)
		outf = open("{}/{}/drilldown.txt".format(indir, tumor_short), "w")
		table = tumor_short+"_simple_somatic"
		outf.write("\n=============================\n")
		outf.write(tumor_short+"\n")
		if tumor_short in cancer_dict: outf.write(cancer_dict[tumor_short]["description"]+"\n")
		outf.write("total donors: {}\n".format(len(get_donors(cursor, table))))
		outf.write("\n")
		gene_groups = parse("{}/{}/{}".format(indir, tumor_short, 'pathways.txt'))
		for group in gene_groups:
			retlines = elaborate(cursor, table, reactome_gene_groups, name2reactome_id,  group)
			if retlines: outf.write("\n".join(retlines)+"\n")
		outf.close()
	cursor.close()
	db.close()

	return


#########################################
if __name__ == '__main__':
	main()
