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


# to further classify tumor subtypes, the following files might be useful:
# specimen*tsv contains tumor histological type code, the ontology can be found here: http://codes.iarc.fr/codegroup/2
# donor*tsv contains donor diagnosis icd10, the ontology can be found here: http://apps.who.int/classifications/icd10/browse/2016/en

import json
import os, subprocess

from config import Config

DEVNULL = open(os.devnull, 'wb')

#########################################
def parse_json(release):

	inf = open ("dirstruct.json", "r")
	jsonstr = inf.read()
	inf.close()

	json_parsed = json.loads(jsonstr)

	# I have a hunch there should be a simpler way
	release_dir = [i for i in json_parsed if i['name']=="/release_%s"%release][0]
	projects = [i for i in release_dir['contents'] if ("Projects" in i['name'])][0]

	of_interest = {}
	for i in projects['contents']:
		project_name = i['name'].split('/')[-1]
		if 'contents' not in i: continue
		has_somatic = False
		names = []
		for filedescr in i['contents']:
			if 'somatic' in filedescr['name']: has_somatic = True
			names.append(filedescr['name'])
		if has_somatic:
			of_interest[project_name] = names

	return of_interest

#########################################
def main():
	icgc_token = os.environ['ICGC_TOKEN']
	base_url = "https://dcc.icgc.org/api/v1/download?fn="

	projects_of_interest = parse_json(Config.icgc_release)
	for project_name, fnms in projects_of_interest.items():
		print(project_name)
		target_dir = "/".join([Config.data_home_local] + project_name.split('-'))
		print(target_dir)
		if not os.path.exists(target_dir): os.makedirs(target_dir)
		for fnm_full_path in fnms:
			if not 'somatic' in fnm_full_path: continue
			# not sure how to download the listing of controlled files
			# I am going with the assumption that for every 'open' file
			# there is a 'controlled' version
			fnm_full_path = fnm_full_path.replace('open','controlled')
			fnm = fnm_full_path.split('/')[-1]
			print(("\t", fnm))
			if os.path.exists("/".join([target_dir, fnm])): continue
			if fnm_full_path[0]=='/': fnm_full_path=fnm_full_path[1:]
			cmd = "curl -L '{}/{}' -o {}/{} --header 'authorization: Bearer {}' ".\
				format(base_url, fnm_full_path, target_dir, fnm, icgc_token)
			subprocess.call(["bash", "-c", cmd], stdout=DEVNULL, stderr=DEVNULL)

#########################################
if __name__ == '__main__':
	main()
