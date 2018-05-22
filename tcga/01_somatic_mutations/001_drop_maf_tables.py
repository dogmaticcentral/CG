#!/usr/bin/python

#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2016 Ivana Mihalek.
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
import subprocess

from tcga_utils.mysql import *


#########################################
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	tcga_home = "/storage/databases/tcga"
	cmd = "find %s -name Somatic_Mutations" % tcga_home
	cancer_types = [dirnm.split("/")[4] for dirnm in subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE).stdout.readlines()]

	# I don't want to start this by mistake - remove comment and put in db names as needed
	print "please comment out if you are sure you want to delete database tables"
	exit(1)

	db_name = "tcga"
	switch_to_db(cursor,db_name)

	for cancer_type in cancer_types:
		# check somatic mutations table exist
		#for table in ( '%s_somatic_mutations'%cancer_type, '%s_mutations_meta'%cancer_type, '%s_conflict_mutations'%cancer_type):
		for table in ( '%s_somatic_mutations'%cancer_type, '%s_conflict_mutations'%cancer_type):
			if (check_table_exists (cursor, db_name, table)):
				print table, " found in ", db_name
				qry = "drop table %s "  % table
				rows = search_db(cursor, qry)
			else:
				print table, " not found in ", db_name

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
