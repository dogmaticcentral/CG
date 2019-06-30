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
from time import time
import re

# use R functions - rpy2
from rpy2.robjects.packages import importr
# in particular, p-value under hypergeometric distribution
phyper = importr('stats').phyper

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


#########################################
def silent_nonsilent_retrieve(cursor, gene):
	transcript_id = approved_symbol2ensembl_canonical_transcript(cursor, gene)
	if not transcript_id: return -5, -5
	qry = "select silent, nonsilent from ensembl_coding_seqs where transcript_id='%s' " % transcript_id
	ret = error_intolerant_search(cursor, qry)
	if ret:
		return ret[0]
	else:
		return -6, -6

#
#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./icgc_utils/kernprof.py -l 75_....py
# followed by
# python3 -m line_profiler 75_....py.lprof
# see here https://github.com/rkern/line_profiler#line-profiler
# the reason I am using local kernprof.py is that I don't know where pip
# installed its version (if anywhere)
#@profile
####################################################
def elaborate (cursor, table, reactome_gene_groups, name2reactome_id, group, mut_type_characterization):
	[name, number_of_genes, cdna_length, donors_affected, expected_donors, stdev, zscore] = group
	if abs(zscore)<5: return None
	reactome_id = name2reactome_id[name]
	retlines = []
	retlines.append("{}  {}".format(name, zscore))
	retlines.append("\t number of genes:  %3d   combined cdna length: %d " % (number_of_genes, cdna_length) )
	retlines.append("\t number of donors affected:  %3d " %  donors_affected)
	retlines.append("\t expected donors affected:  %6.2f      stdev:  %6.2f    z:  %6.2f " % (expected_donors, stdev, zscore))
	if donors_affected==0: return

	retlines.append("\tall genes in the Reactome group %d " % len(reactome_gene_groups[reactome_id]))
	retlines.append("\t"+str(sorted(reactome_gene_groups[reactome_id])))
	gene_string = ",".join([quotify(g) for g in reactome_gene_groups[reactome_id]])
	qry  = "select  gene_symbol,  count(distinct(icgc_donor_id)) from %s " % table
	qry += "where pathogenicity_estimate=1 and reliability_estimate=1 "
	qry += "and gene_symbol in (%s) group by gene_symbol" % gene_string
	# Under Python 3.6, the built-in dict does track insertion order,
	# although this behavior is a side-effect of an implementation change and should not be relied on.
	genes_mutated = dict(sorted(error_intolerant_search(cursor, qry), key= lambda r: r[1], reverse=True))
	retlines.append("\tactually affected %d " % len(genes_mutated))
	retlines.append("\tdonors per affected gene: "+str(genes_mutated))
	for gene in genes_mutated.keys():
		if gene in mut_type_characterization:
			retlines.extend(mut_type_characterization[gene])
			continue
		temp_storage = []
		silent_possible, nonsilent_possible =  silent_nonsilent_retrieve(cursor, gene)
		temp_storage.append("\t %s all possible:  silent %d    nonsilent %d    silent/nonsilent %.3f    silent/(silent+nonsilent) %.3f " %
						(gene, silent_possible, nonsilent_possible, silent_possible/nonsilent_possible,
						 silent_possible/(silent_possible+nonsilent_possible) ))
		expected_silent_to_nonsilent = silent_possible/nonsilent_possible
		chromosome = find_chromosome(cursor, gene)
		qry  = "select m.mutation_type, count(*) c from  %s s, mutations_chrom_%s m where s.gene_symbol = '%s' " % (table, chromosome, gene)
		qry += "and s.icgc_mutation_id=m.icgc_mutation_id  "
		qry += "and s.pathogenicity_estimate=1 and s.reliability_estimate=1  "
		qry += "group by m.mutation_type"
		ret = error_intolerant_search(cursor, qry)
		if not ret or len(ret)==0: continue
		type_count = dict(ret)
		for muttype, ct in type_count.items():
			temp_storage.append("\t\t {}: {}".format(muttype, ct))
			if muttype != 'single': continue
			qry = "select m.consequence, count(*) c from  %s s, mutations_chrom_%s m where s.gene_symbol = '%s' " % (table, chromosome, gene)
			qry += "and s.icgc_mutation_id=m.icgc_mutation_id  and m.mutation_type='single'  "
			# we want silent mutations here, that are already marked with pathg 0 (thus drop pathg=1 filter)
			qry += "and s.reliability_estimate=1  "
			qry += "group by m.consequence"
			ret = error_intolerant_search(cursor, qry)
			if not ret or len(ret)==0: continue
			conseq_ct = dict(ret)
			nonsilent = conseq_ct.get('missense',0) + conseq_ct.get('stop_gained',0) + conseq_ct.get('stop_lost',0)
			silent = conseq_ct.get('synonymous',0)
			if silent+nonsilent==0: continue
			for conseq, ct in conseq_ct.items():
				temp_storage.append("\t\t\t\t {}: {}".format(conseq, ct))
			silent_to_nonsilent  = float(silent)/nonsilent if nonsilent>0 else 1.0
			# pvalue under hypergeometric distribution (R)
			# https://stackoverflow.com/questions/8382806/hypergeometric-test-phyper
			# phyper(success in sample, success in population, failure in population, sample size)
			# that is
			# phyper(silent, silent_possible, nonsilent_possible, silent+nonsilent)
			phyper_ret = phyper(silent, silent_possible, nonsilent_possible, silent+nonsilent)
			pval = phyper_ret[0] if phyper_ret else 1.0
			temp_storage.append("\t\t silent/nonsilent = %.3f    expected = %.3f    pval = %.0e" %
			                    (silent_to_nonsilent, expected_silent_to_nonsilent, pval))
		if len(temp_storage)>0:
			mut_type_characterization[gene] = temp_storage
			retlines.extend(mut_type_characterization[gene])
	retlines.append("")
	return retlines


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
	#for tumor_short in ['AML']:
		print(tumor_short)
		time0 = time()
		outf = open("{}/{}/drilldown.txt".format(indir, tumor_short), "w")
		table = tumor_short+"_simple_somatic"
		outf.write("\n=============================\n")
		outf.write(tumor_short+"\n")
		if tumor_short in cancer_dict: outf.write(cancer_dict[tumor_short]["description"]+"\n")
		outf.write("total donors: {}\n".format(len(get_donors(cursor, table))))
		outf.write("\n")
		gene_groups = parse("{}/{}/{}".format(indir, tumor_short, 'pathways.txt'))
		mut_type_characterization = {} # store results for individual genes
		for group in gene_groups:
			retlines = elaborate(cursor, table, reactome_gene_groups, name2reactome_id,  group, mut_type_characterization)
			if retlines: outf.write("\n".join(retlines)+"\n")
		outf.close()
		print("\t done in  %.2f mins" % (float(time()-time0)/60))
	cursor.close()
	db.close()

	return


#########################################
if __name__ == '__main__':
	main()
