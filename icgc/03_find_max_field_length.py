#! /usr/bin/python3

# Note: this is really slow, but assumed to be run only once - during the installation

# for the  ICGC v27 the search for the longest entry  comes up with the following

# alleles 200
# genotypes 401
# gene_names 15
# aa_mutation 59
# cds_mutation 12


#the longest genotype I see in this set is 401
#I assume that the allele legnth cutoff is 200
# I'll make the allele columns varchar210 and genotype 430

import os
from config import Config

#########################################
def get_simple_somatic_tsv_files(data_home_local):
	tsv_files = []
	for root, dirs, files in os.walk(data_home_local):
		for file in files:
			if file.endswith(".tsv") and 'simple_somatic' in file:
				tsv_files.append(os.path.join(root, file))
	return tsv_files


#########################################
def main():

	tsv_files = get_simple_somatic_tsv_files(Config.data_home_local)

	field_groups = {
		'alleles':      ['reference_genome_allele', 'mutated_from_allele', 'mutated_to_allele'],
		'genotypes':    ['control_genotype', 'tumour_genotype', 'control_genotype'],
		'gene_names':   ['gene_affected', 'transcript_affected'],
		'aa_mutation':  ['aa_mutation'],
		'cds_mutation': ['cds_mutation']
	}

	max_allele_length = {}
	for field_group in field_groups.keys():
		max_allele_length[field_group] = 0

	for tf in tsv_files:
		print(tf)
		infile = open(tf,'r')
		headers = None
		for line in infile:
			if not headers:
				headers = line.rstrip('\n').split('\t')
			else:
				fields = line.rstrip('\n').split('\t')
				field_named = dict(list(zip(headers, fields)))
				for field_group, field_names in field_groups.items():
					for name in field_names:
						if len(field_named[name])>max_allele_length[field_group]:
							max_allele_length[field_group] = len(field_named[name])
		infile.close()

	for field_group in field_groups.keys():
		print(field_group, max_allele_length[field_group])

#########################################
if __name__ == '__main__':
	main()
