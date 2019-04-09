#!/usr/bin/python3
#
#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2019 Ivana Mihalek.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
#
# Contact: ivana.mihalek@gmail.com
#


import requests
# Set security warning to always go off by default.
import json
import sys, os


########################################
def main():

	# working with legacy from TCGA
	# the list of maf files downloaded from
	# https://portal.gdc.cancer.gov/legacy-archive/search/f?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_category%22,%22value%22:%5B%22Simple%20nucleotide%20variation%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_format%22,%22value%22:%5B%22MAF%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D%5D%7D
	# simple somatic files, maf format, open access
	json_data = json.loads(open("/data/tcga/files.2018-05-20.json").read())
	# the files proper should be in /data/tcs/dumpspace, downloaded by the download client provided by gdc
	# a problem: some files are revised versions of older files - htat info is not available

	files_per_cancer_type = {}

	for filedata in json_data:
		#print filedata['file_name'], filedata['file_id']
		project_ids = set([case['project']['project_id'] for case in filedata['cases']])
		#print "\t", project_ids
		if len(project_ids)!=1:
			print "unexpected length of project ids",  len(project_ids), filedata['file_name']
		cancer_type = list(project_ids)[0].split("-")[1]
		if not files_per_cancer_type.has_key(cancer_type):
			files_per_cancer_type[cancer_type] = []
		files_per_cancer_type[cancer_type].append(filedata['file_id'])
		#file_id[filedata['file_name']] = filedata['file_id']


	for cancer_type, files in files_per_cancer_type.iteritems():

		target_path = "/data/tcga/%s/Somatic_Mutations" % cancer_type

		print cancer_type
		print target_path
		if not os.path.exists(target_path):
			print target_path, "not found - creating"
			os.makedirs(target_path)

		name_seen = {}
		for file_id in files_per_cancer_type[cancer_type]:

			source_path = "/data/tcga/dumpspace/{}".format(file_id)
			if not os.path.exists(source_path):
				print "\t", source_path, "not found"
				continue

			for path, dirs, files in os.walk(source_path):
				mafs = [fn for fn in files if fn[-3:]=="maf"]
				#if len(mafs)>0:  print path, "   ".join(mafs)
				for maf in mafs:
					if not name_seen.has_key(maf): name_seen[maf] = 0
					name_seen[maf] += 1
					if name_seen[maf]>1:
						new_name = "{}.{}.maf".format(maf[:-4],name_seen[maf])
					else:
						new_name = maf
					src  = "%s/%s" % (path,maf)
					dest = "%s/%s" % (target_path, new_name)
					print "move   %s   to  %s " % (src, dest)
					os.rename(src,dest)


			#duplicated_names = [name for name in name_seen.keys() if name_seen[name]>1]
			#if len(duplicated_names)>0: print  " ... ", "  ".join(duplicated_names)

		#exit()

# unzip: find /mnt/databases/TCGA -name "*maf.gz" | xargs gunzip

#########################################
if __name__ == '__main__':
	main()

