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


import time

from config import Config
from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *


def delete_normal(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, 'icgc')
	for table in tables:
		print("==================")
		print(table)
		tumor = table.split("_")[0]
		qry  = "select icgc_specimen_id from %s_specimen " % tumor
		qry += "where specimen_type like 'Normal%'"
		ret = search_db(cursor, qry, verbose=False)
		if not ret: continue
		normal_specimen_ids = [r[0] for r in ret]
		for icgc_specimen_id in normal_specimen_ids:
			qry  = "delete from %s " % table
			qry += "where  icgc_specimen_id='%s' " % icgc_specimen_id
			search_db(cursor, qry, verbose=False)

	cursor.close()
	db.close()


#########################################
#########################################

def main():

	print("disabled")
	exit()

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	table_size = get_table_size(cursor, 'icgc', tables)
	cursor.close()
	db.close()

	tables_sorted = sorted(tables, key=lambda t: table_size[t], reverse=True)
	half = int(len(tables_sorted)/2)
	tables_mirrored  = tables_sorted[0:half] + list(reversed(tables_sorted[half:]))
	number_of_chunks = int(half/2)

	parallelize(number_of_chunks, delete_normal, tables_mirrored, [], round_robin=True)



#########################################
if __name__ == '__main__':
	main()
