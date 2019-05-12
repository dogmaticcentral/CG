#!/usr/bin/python
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
from config import Config
from icgc_utils.mysql   import  *


#########################################
#########################################
def main():

	print("disabled")
	return

	cancer_types = []
	for name in os.listdir(Config.data_home_local):
		if os.path.isdir("/".join([Config.data_home_local,name])): cancer_types.append(name)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	db_name = "icgc"
	switch_to_db(cursor, db_name)
	for ct in cancer_types: # note we are filling the temp table
		for filetype in ["simple_somatic_temp", "donor", "specimen"]:
			table = "{}_{}".format(ct, filetype)
			qry = "load data local infile 'tsvs/%s.tsv' into table %s" % (table,table)
			search_db(cursor,qry,verbose=True)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()

