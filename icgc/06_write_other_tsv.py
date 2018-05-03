#! /usr/bin/python

# just spit out the select columns ans slurp them into database

import os, subprocess
from subprocess import PIPE

# import the prorudeced fiels using mysqlimport
# mysqlimport db_name table_name.ext
# by convention, mysql takes that the input file name is the table_name.ext
# mysqlimport will strip the extension by itself
# the default separator is \t, so the .ext will typically be something like tsv or csv
# => provided there are no unrelated tsv's in the workdir, something like
# for i in `ls *.tsv`; do mysqlimport -u name
# should do
# fi you get the --secure=file-priv crap try


#########################################
def get_tsv_files(data_home_local, type):
	tsv_files = []
	for root, dirs, files in os.walk(data_home_local):
		for file in files:
			if file.endswith(".tsv") and type+"." in file:
				tsv_files.append(os.path.join(root, file))
	return tsv_files

#########################################
def appendopen(original_tsv_file, filetype):
	fields = original_tsv_file.split("/")
	cancer_type  = fields[3]
	outname = "{}_{}.tsv".format(cancer_type, filetype)
	if os.path.exists(outname):
		last_id = int(subprocess.Popen(["bash", "-c", "tail -n1 %s"%outname], stdin=PIPE, stdout=PIPE).communicate()[0].split("\t")[0])
	else:
		last_id = 0
	print "writing to", outname, "prev id:", last_id

	return open(outname,'a'), last_id

#########################################
def main():
	data_home_local = "/data/icgc"

	names_string = {"donor":"icgc_donor_id,submitted_donor_id,donor_sex,donor_diagnosis_icd10",
	                "specimen":"icgc_specimen_id,icgc_donor_id,specimen_type,tumour_histological_type"}

	for filetype in ["donor", "specimen"]:
		tsv_files = get_tsv_files(data_home_local, filetype)
		names = names_string[filetype].split(",")

		outfiles = []
		for tf in tsv_files:
			print tf
			infile  = open(tf,'r')
			outfile, last_id = appendopen(tf,filetype)
			if not outfile in outfiles: outfiles.append(outfile)
			headers = None
			id = last_id
			for line in infile:
				if not headers:
					headers = line.rstrip('\n').split('\t')
				else:
					fields = line.rstrip('\n').split('\t')
					field_named = dict(zip(headers, fields))
					id += 1
					new_fields = [str(id)]
					for name in names:
						new_fields.append(field_named[name])
					outfile.write("\t".join(new_fields) + "\n")
			infile.close()
			outfile.close()


#########################################
if __name__ == '__main__':
	main()
