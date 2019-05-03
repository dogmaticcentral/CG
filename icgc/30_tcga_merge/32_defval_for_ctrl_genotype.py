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
from config import Config
from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *

#########################################
#########################################
def main():

	#print("disabled")
	#exit()

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,"icgc")
	for table in tables:
		print("\n====================")
		print(table)
		# getting rid of "default null"
		#  perhaps not needed if not set in 10_local_db_loading/06_make_tables.py
		qry = "ALTER TABLE %s modify column  control_genotype  VARCHAR (430)" % table
		search_db(cursor, qry, verbose=True)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
