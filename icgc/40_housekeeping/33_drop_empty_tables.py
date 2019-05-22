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


from icgc_utils.common_queries   import  *
from config import Config


#########################################
#########################################
# In the cleaning process we might have ended upe with some empty tables.
# In the ICGC v27, for example, all NACA tables consisted of normal donors
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")
	for table in tables:
		number_of_rows = hard_landing_search(cursor, "select count(*) from %s"%table)[0][0]
		if number_of_rows>0: continue
		print(table, number_of_rows)
		error_intolerant_search(cursor, "drop table %s"%table)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
