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
from icgc_utils.processes import *


def cleanup_duplicate_specimens(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for somatic_table in tables:
		ret = search_db(cursor, "select distinct icgc_donor_id from %s"%somatic_table)
		if not ret: continue
		icgc_donor_ids = [r[0] for r in ret]
		problematic = []
		specimen_ids = {}
		for icgc_donor_id in icgc_donor_ids:
			qry  = "select distinct(icgc_specimen_id) "
			qry += "from %s where icgc_donor_id='%s' " % (somatic_table,icgc_donor_id)
			ret  = search_db(cursor,qry)
			if ret and len(ret)>1:
				problematic.append(icgc_donor_id)
				specimen_ids[icgc_donor_id] = [r[0] for r in ret]
		if len(problematic)==0:
			print("%s has no duplicates in the variant table" % somatic_table)
			continue

		print("%s has %d donor ids with duplicate specimen ids in the variant table" % (somatic_table,len(problematic)))

		# most of these are innocuous, with normal sample not appearing in the variants table
		# this, however is not always the case
		resolve_duplicate_specimens(cursor, somatic_table, specimen_ids)

	cursor.close()
	db.close()


	return


#########################################
#########################################
def main():

	print("disabled - this script deletes certain rows ") # comment out to run
	exit(1)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which  somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]

	cursor.close()
	db.close()

	number_of_chunks = 1
	parallelize(number_of_chunks, cleanup_duplicate_specimens, tables, [])

	return



#########################################
if __name__ == '__main__':
	main()
