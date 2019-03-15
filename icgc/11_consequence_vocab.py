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

import MySQLdb

from config import Config
from icgc_utils.mysql   import  *

#########################################
#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,'icgc')
	cons_vocab = set()
	for mutations_table in tables:
		print(mutations_table)
		qry = "select distinct consequence_type from %s" % mutations_table
		csq = [ret[0] for ret in search_db(cursor,qry)]
		cons_vocab = cons_vocab.union(csq)

	for kwd in cons_vocab:
		print(kwd)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()



