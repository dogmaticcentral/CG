#! /usr/bin/python


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
				for name in ['sequencing_strategy', 'quality_score', 'probability', 'total_read_count']:
						print name +":", field_named[name],
				print

		infile.close()

#########################################
if __name__ == '__main__':
	main()
