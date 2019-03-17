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
from icgc_utils.icgc   import  *


#################
def make_location_table(cursor, table):
	# location:  chromosome_position  gene_relative[semicolon separated list of the format ENSG:relative_location]
	# transcript_relative[semicolon separated list of the format ENST:relative_location]

	qry = "drop table " + table
	search_db(cursor, qry, verbose=True)
	qry = ""
	qry += "  CREATE TABLE  %s (" % table
	qry += "	 position INT  NOT NULL, "
	qry += "     gene_relative TEXT, "
	qry += "     transcript_relative TEXT, "

	qry += "	 PRIMARY KEY (position) " # automatically indexed
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)
	return


#########################################
#########################################
def main():

	print("disabled")
	exit()

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")

	# indices on simple somatic temp
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like 'locations_%'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	for table in tables:
		print(table)
		make_location_table(cursor, table)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()

