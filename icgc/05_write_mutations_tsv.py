#! /usr/bin/python

# just spit out the select columns ans slurp them into database

import os, subprocess
from subprocess import PIPE

# import the produced fields using mysqlimport
# mysqlimport db_name table_name.ext
# by convention, mysql takes that the input file name is the table_name
# mysqlimport will strip the extension by itself
# the default separator is \t, so the .ext will typically be something like tsv or csv
# => provided there are no unrelated tsv's in the workdir, something like
# for i in `ls *.tsv`; do mysqlimport -u name
# should do the trick
# if you get the --secure=file-priv crap try https://stackoverflow.com/questions/32737478/how-should-i-tackle-secure-file-priv-in-mysql


#########################################
from config import Config


def get_simple_somatic_tsv_files(data_home_local):
	tsv_files = []
	for root, dirs, files in os.walk(data_home_local):
		for file in files:
			if file.endswith(".tsv") and 'simple_somatic' in file:
				tsv_files.append(os.path.join(root, file))
	return tsv_files

#########################################
def appendopen(original_tsv_file):
	if not os.path.exists("tsvs"): os.mkdir("tsvs")
	# the first thing after the storage path should be the cancer name
	cancer_type = original_tsv_file[len(Config.data_home_local)+1:].split("/")[0]
	outname = "tsvs/"+cancer_type+"_simple_somatic_temp.tsv"
	if os.path.exists(outname):
		last_id = int(subprocess.Popen(["bash", "-c", "tail -n1 %s"%outname], stdin=PIPE, stdout=PIPE).communicate()[0].split("\t")[0])
	else:
		last_id = 0
	print("writing to", outname, "prev id:", last_id)

	return open(outname,'a'), last_id

#########################################
def main():

	tsv_files = get_simple_somatic_tsv_files(Config.data_home_local)

	# get this by
	# head -n1 simple_somatic_mutation.controlled.ALL-US.tsv | sed 's/\t/,/g'
	# (for example)
	names = "icgc_mutation_id,icgc_donor_id,project_code,icgc_specimen_id," \
			"icgc_sample_id,matched_icgc_sample_id,submitted_sample_id," \
			"submitted_matched_sample_id,chromosome,chromosome_start,chromosome_end," \
			"chromosome_strand,assembly_version,mutation_type,reference_genome_allele,control_genotype," \
			"tumour_genotype,expressed_allele,mutated_from_allele,mutated_to_allele,quality_score,probability," \
			"total_read_count,mutant_allele_read_count,verification_status,verification_platform," \
			"biological_validation_status,biological_validation_platform,consequence_type,aa_mutation," \
			"cds_mutation,gene_affected,transcript_affected,gene_build_version,platform,experimental_protocol," \
			"sequencing_strategy,base_calling_algorithm,alignment_algorithm,variation_calling_algorithm," \
			"other_analysis_algorithm,seq_coverage,raw_data_repository," \
			"raw_data_accession,initial_data_release_date".split(",")

	outfiles = []

	for tf in tsv_files:
		print(tf)
		infile  = open(tf,'r')
		outfile, last_id = appendopen(tf)
		if not outfile in outfiles: outfiles.append(outfile)
		headers = None
		id = last_id

		for line in infile:
			if not headers:
				headers = line.rstrip('\n').split('\t')
			else:
				fields = line.rstrip('\n').split('\t')
				field_named = dict(list(zip(headers, fields)))
				id += 1
				new_fields = [str(id)]
				# some fixing
				field_named['mutation_type']    = field_named['mutation_type'].split(" ")[0]
				field_named['consequence_type'] = field_named['consequence_type'].replace("_variant","").replace("_gene","")
				for name in ["total_read_count","mutant_allele_read_count"]:
					tmp = field_named[name].replace(" ","")
					if len(tmp)==0: tmp = "\N"
					field_named[name] = tmp
				for name in names:
					new_fields.append(field_named[name])
				outfile.write("\t".join(new_fields) + "\n")

		infile.close()
		outfile.close()


#########################################
if __name__ == '__main__':
	main()
