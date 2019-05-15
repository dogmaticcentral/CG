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
import time
from icgc_utils.processes   import  *

def cleanup(tables, other_args):
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for somatic_table in tables:
		print("\n========================================")
		print(somatic_table)
		# creating indices at this point makes sense only if we expect to run this repeatedly
		# this won't do anything if the index exists already
		create_index (cursor, 'icgc', 'spec_idx', somatic_table, ['icgc_specimen_id'])
		print( " created spec_idx on %s  " % somatic_table )
		qry  = "select icgc_specimen_id, count(distinct(icgc_sample_id)) as c "
		qry += "from %s group by icgc_specimen_id having c>1 " % somatic_table
		ret = search_db(cursor,qry)
		if not ret:
			print("%s has no duplicate samples in the variant table"% somatic_table)
			continue

		problematic = [r[0] for r in ret]
		print("%s has %d specimen ids with duplicate sample ids in the variant table" % (somatic_table,len(problematic)))
		continue
		# most of these are innocuous, with normal sample not appearing in the variants table
		# this, however is not always the case
		for icgc_specimen_id in problematic:
			# however, the mutations are not necessarily  duplicated ...
			qry  = "select icgc_mutation_id, count(*) as c "
			qry += "from %s where  icgc_specimen_id='%s' " % (somatic_table,icgc_specimen_id)
			qry += "group by icgc_mutation_id having c>1 "
			ret2 = search_db(cursor,qry)
			if not ret2 or len(ret2)==0:
				print("\t %s no duplicate mutation ids" %icgc_specimen_id)
			else:
				duplicates = dict(ret2)
				for icgc_mutation_id, ct in duplicates.items():
					qry  = "select * from %s " % somatic_table
					qry += "where icgc_specimen_id='%s' and  icgc_mutation_id='%s'" % (icgc_specimen_id, icgc_mutation_id)
					ret3 = search_db(cursor,qry)
					resolve_duplicate_mutations(cursor, somatic_table, ret3, verbose=False)

	cursor.close()
	db.close()


#########################################
#########################################
def main():

	#print("disabled ") # comment out to run
	#exit(1)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which  somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	table_size = get_table_size(cursor, 'icgc', tables)
	cursor.close()
	db.close()

	tables_sorted = sorted(tables, key=lambda t: table_size[t], reverse=True)
	half = int(len(tables_sorted)/2)
	tables_mirrored  = tables_sorted[0:half] + list(reversed(tables_sorted[half:]))
	number_of_chunks = half

	number_of_chunks = 1
	tables_mirrored = ['THCA_simple_somatic']

	parallelize(number_of_chunks, cleanup, tables_mirrored, [], round_robin=True)


	return



#########################################
if __name__ == '__main__':
	main()
