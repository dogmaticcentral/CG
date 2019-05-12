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

# just spit out the select columns ans slurp them into database

import os, subprocess
from subprocess import PIPE
from icgc_utils.mysql import *
from icgc_utils.processes import *

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


def get_simple_somatic_tsv_files(data_home_local, cancer_types):
	tsv_files = []
	outstr = ""
	for ct in cancer_types:
		outstr += ct + "\n"
		for root, dirs, files in os.walk("{}/{}".format(data_home_local,ct)):
			for file in files:
				if file.endswith(".tsv") and 'simple_somatic' in file:
					outstr += "\t" + os.path.join(root, file)  + "\n"
					tsv_files.append(os.path.join(root, file))
		print("{} writing:\n{}".format(os.getpid(), outstr))
	return tsv_files

#########################################
def appendopen(original_tsv_file):
	# the first thing after the storage path should be the cancer name
	cancer_type = original_tsv_file[len(Config.data_home_local)+1:].split("/")[0]
	outname = "tsvs/"+cancer_type+"_simple_somatic_temp.tsv"
	if os.path.exists(outname):
		last_id = int(subprocess.Popen(["bash", "-c", "tail -n1 %s"%outname], stdin=PIPE, stdout=PIPE).communicate()[0].decode('utf_8').split("\t")[0])
	else:
		last_id = 0
	print("writing to", outname, "prev id:", last_id)

	return open(outname,'a'), last_id

#########################################
def write_tsvs(cancer_types, other_args):


	tsv_files = get_simple_somatic_tsv_files(Config.data_home_local, cancer_types)

	# the full set of the header fields can be found by (for example)
	# head -n1 simple_somatic_mutation.controlled.ALL-US.tsv | sed 's/\t/,/g'
	# this correposnds to the subset we want to store
	names = "icgc_mutation_id,icgc_donor_id,icgc_specimen_id,icgc_sample_id,submitted_sample_id,chromosome,chromosome_start," \
	"chromosome_end,chromosome_strand,assembly_version,mutation_type,reference_genome_allele,control_genotype,tumour_genotype,"\
	"mutated_from_allele,mutated_to_allele,consequence_type,aa_mutation,cds_mutation," \
	"gene_affected,transcript_affected,total_read_count,mutant_allele_read_count".split(",")

	for tf in tsv_files:
		print(tf)
		infile  = open(tf,'r')
		outfile, last_id = appendopen(tf)
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
					if len(tmp)==0: tmp = "\\N"
					field_named[name] = tmp
				for name in names:
					new_fields.append(field_named[name])
				outfile.write("\t".join(new_fields) + "\n")

		infile.close()
		outfile.close()

#########################################
def main():

	print("disabled")
	exit()

	if not os.path.exists("tsvs"): os.mkdir("tsvs")

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")

	# indices on simple somatic temp
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	cancer_types = [field[0].split("_")[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	number_of_chunks = 10  # myISAM does not deadlock
	parallelize(number_of_chunks, write_tsvs, cancer_types, [])



#########################################
if __name__ == '__main__':
	main()
