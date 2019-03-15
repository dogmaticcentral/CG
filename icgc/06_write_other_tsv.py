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
# note: in v 27 in PCSI nad in COCA CN in handful of places there  are inconsistencies
# icd10 diagnosis annotation (use of quotes, some kind of extended annotation)
# affected ids: COCA-CN-LPG,  COCA-CN-NYG, COCA-CN-WQ, COCA-CN-CJ, COCA-CN-JKQ, PCSI_0506
#
import os, subprocess
from subprocess import PIPE


#########################################
from config import Config


def get_tsv_files(data_home_local, type):
	tsv_files = []
	for root, dirs, files in os.walk(data_home_local):
		for file in files:
			if file.endswith(".tsv") and type+"." in file:
				tsv_files.append(os.path.join(root, file))
	return tsv_files

#########################################
def appendopen(original_tsv_file, filetype):
	if not os.path.exists("tsvs"): os.mkdir("tsvs")
	# the first thing after the storage path should be the cancer name
	cancer_type = original_tsv_file[len(Config.data_home_local)+1:].split("/")[0]
	outname = "tsvs/{}_{}.tsv".format(cancer_type, filetype)
	if os.path.exists(outname):
		last_id = int(subprocess.Popen(["bash", "-c", "tail -n1 %s"%outname], stdin=PIPE, stdout=PIPE).communicate()[0].decode('utf_8').split("\t")[0])
	else:
		last_id = 0
	print("writing to", outname, "prev id:", last_id)

	return open(outname,'a'), last_id

#########################################
def main():

	names_string = {"donor":"icgc_donor_id,submitted_donor_id,donor_sex,donor_diagnosis_icd10",
	                "specimen":"icgc_specimen_id,icgc_donor_id,specimen_type,tumour_histological_type"}

	for filetype in ["donor", "specimen"]:
		tsv_files = get_tsv_files(Config.data_home_local, filetype)
		names = names_string[filetype].split(",")

		outfiles = []
		for tf in tsv_files:
			print(tf)
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
					field_named = dict(list(zip(headers, fields)))
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
