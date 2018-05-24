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

import os.path
import re, commands, time
from tcga_utils.processes import  *
from tcga_utils.utils import  *
import time


##################################################################################
def main():

	db_dir = '/storage/databases/tcga'
	if not os.path.isdir(db_dir):
		print "directory " + db_dir + " not found"
		exit(1)  # TCGA db dir not found

	db     = connect_to_mysql()
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='tcga' and table_name like '%_somatic_mutations'"
	tables = [field[0] for field in search_db(cursor,qry)]

	switch_to_db(cursor,"tcga")
	for table in tables:
		qry = "create index location_idx on %s (chromosome,start_position,end_position)" % table
		search_db(cursor,qry,verbose=True)
	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
