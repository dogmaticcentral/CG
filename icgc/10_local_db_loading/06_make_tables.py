#!/usr/bin/python3
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

from icgc_utils.mysql import  *
from icgc_utils.icgc  import  *
from config import Config

#########################################
#########################################
def main():

	print("disabled bcs it drops and re-creates tables - comment out to run")
	return

	homedir = Config.data_home_local
	cancer_types = []
	for name in os.listdir(homedir): # this is not walking recursively - it's just this dir
		if os.path.isdir("/".join([homedir,name])): cancer_types.append(name)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	db_name =  "icgc"
	for ct in cancer_types:
		# note 'temp' here - we'll reroganize when we get to TCGA
		mutations_table = ct + "_simple_somatic_temp"
		make_temp_somatic_muts_table(cursor, db_name, mutations_table)
		donors_table = ct + "_donor"
		make_donors_table(cursor, db_name, donors_table)
		specimen_table = ct + "_specimen"
		make_specimen_table(cursor, db_name, specimen_table)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()

