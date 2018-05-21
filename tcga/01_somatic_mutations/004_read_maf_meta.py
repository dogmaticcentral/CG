#!/usr/bin/python -u
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
# store the meta info about the maf files: name and the reference genome,  for now

import os
from old_tcga_tools.tcga_utils.mysql import *
from old_tcga_tools.tcga_utils.utils import *
from random import random
import urllib2
from HTMLParser import HTMLParser
from bs4 import BeautifulSoup

##################################################################################
##################################################################################
def main():
	db = connect_to_mysql()
	cursor = db.cursor()

	qry = "select table_name from information_schema.tables "
	qry += "where table_schema='tcga' and table_name like '%_mutations_meta'"
	tables = [field[0] for field in  search_db(cursor,qry)]

	switch_to_db(cursor, "tcga")

	for table in tables:

		print table
		qry = "select * from %s" % table
		rows = search_db(cursor,qry)
		if not rows:
			print "\t no meta info found"
			continue
		for row in rows:
			[meta_id, file_name, quality_check, assembly, diagnostics] = row
			#if diagnostics and "tumor alleles identical" in diagnostics:
			if True:
				print "\t %4d  %50s   " % (meta_id, file_name)
				print "\t\t quality check: %6s" % (quality_check)
				print "\t\t assembly: %6s" % (assembly)
				if diagnostics:
					print "\t\t diagnostics:"
					for diag in diagnostics.split(';'):
						print "\t\t\t %6s" % diag.strip()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
