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
# profile decorator is for the use with kernprof (a line profiler):
#  ./icgc_utils/kernprof.py -l 18_cleanup_duplicate_entries.py
# followed by
# python3 -m line_profiler 18_cleanup_duplicate_entries.py.lprof
# see here https://github.com/rkern/line_profiler#line-profiler
# the reason I am using local kernprof.py is that I don't know where pip
# installed its version (if anywhere)
#@profile
def remove_duplicates(table_rows, other_args, verbose=False):

	table  = other_args[0]
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")

	# loop over all duplicate entries
	for line in table_rows:
		[mega_id, ct] = line
		if not mega_id: continue # we do not have one of the IDs (happens with TCGA)
		[icgc_mutation_id, icgc_donor_id, icgc_specimen_id, icgc_sample_id] = mega_id.split("_")

		# check the full length of the entry
		qry  = "select * from %s " % table
		qry += "where icgc_mutation_id = '%s' " % icgc_mutation_id
		qry += "and icgc_donor_id = '%s' " % icgc_donor_id
		qry += "and icgc_specimen_id = '%s' " % icgc_specimen_id
		qry += "and icgc_sample_id = '%s' " % icgc_sample_id
		ret2 = search_db(cursor,qry)

		resolve_duplicate_mutations(cursor, table, ret2, verbose)

	cursor.close()
	db.close()

	print ("\tprocess {} exiting; worked on table {}".format(get_process_id(), table))

	return

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
	number_of_chunks = 10 # myISAM does not deadlock
	for table in tables:
		print("\n====================")
		print("inspecting ", table)
		# column names/headers
		# a hack to get all entries that have all relevant ids identical
		qry = "select concat(icgc_mutation_id,'_', icgc_donor_id,'_',icgc_specimen_id,'_',icgc_sample_id) as mega_id, "
		qry += "count(*) as c from %s  group by mega_id having c>1 " % table
		ret  = search_db(cursor,qry)

		if not ret:
			print("\tno duplicates found in", table)
			continue
		print("\t%s has %d duplicates" % (table, len(ret)))
		#print(ret)

		processes = parallelize(number_of_chunks, remove_duplicates, ret, [table])
		if processes: wait_join(processes)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
