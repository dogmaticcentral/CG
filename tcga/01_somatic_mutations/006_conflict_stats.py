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
import re, commands
from old_tcga_tools.tcga_utils.mysql import  *
from old_tcga_tools.tcga_utils.utils import  get_expected_fields, is_useful, make_named_fields
from time import time


#########################################
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='tcga' and table_name like '%_somatic_mutations'"
	tables = [field[0] for field in search_db(cursor,qry)]

	db_name = "tcga"
	switch_to_db(cursor,db_name)

	conflicts = {}
	entries = {}
	for table in tables:
		qry = "select count(1) from %s " % table
		number_of_entries = search_db(cursor, qry)[0][0]
		entries[table] = number_of_entries
		qry = "select count(1) from %s where conflict is not null" % table
		rows = search_db(cursor, qry)
		conflicts[table] = int(rows[0][0])

	tables_sorted = sorted(conflicts, key=conflicts.__getitem__, reverse=True)
	for table in tables_sorted:
		print table, "entries:", entries[table], "conflicts:", conflicts[table]



	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
