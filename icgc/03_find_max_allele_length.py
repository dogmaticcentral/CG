#! /usr/bin/python


# header names:
#head -n1 /data/icgc/OV/AU/simple_somatic_mutation.controlled.OV-AU.tsv \
# | sed 's/\t/\n/g' | awk 'ct +=1 {print ct, " ", $1}'

#the longest genotype I see in this set (aPr 2018) is 401
#I assume that the allele legnth cutoff is 200
# I'll make the allele columns varchar210 and genotype 430

import os, subprocess

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
	data_home_local = "/data/icgc"
	tsv_files = get_simple_somatic_tsv_files(data_home_local)

	max_allele_length = 0
	for tf in tsv_files:
		print tf
		infile = open(tf,'r')
		headers = None
		for line in infile:
			if not headers:
				headers = line.rstrip('\n').split('\t')
			else:
				fields = line.rstrip('\n').split('\t')
				field_named = dict(zip(headers, fields))
				#for name in ['reference_genome_allele', 'control_genotype', 'tumour_genotype',
				#             'control_genotype', 'mutated_from_allele', 'mutated_to_allele']:
				#for name in ['gene_affected', 'transcript_affected']:
				#for name in ['aa_mutation']:
				for name in ['cds_mutation']:
					if len(field_named[name])>max_allele_length:
						max_allele_length = len(field_named[name])
						print name, max_allele_length, field_named[name]

		infile.close()

#########################################
if __name__ == '__main__':
	main()
