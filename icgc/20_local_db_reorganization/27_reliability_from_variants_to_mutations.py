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

# add info about the reliability to mutations_chrom_* tables
# for faster search

import time

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *
from config import Config
from random import shuffle

#########################################
def update_reliability_column_in_mutation_tables(cursor, table):
	chromosomes = [str(i) for i in range(1,23)] + ["X","Y"]
	for chromosome in chromosomes:
		qry  = "select distinct(icgc_mutation_id)  from %s " % table
		qry += "where chromosome='%s' and reliability_estimate=1" % chromosome
		ret = search_db(cursor,qry, verbose=False)
		if not ret: continue

		reliable_muts = [r[0] for r in ret]
		for batch in range(len(reliable_muts),100):
			batchstring =  ",".join(["'%s'"%m for m in reliable_muts[batch:batch+100]])
			qry  = "update %s " % table
			qry += "set reliability_estimate=1 where icgc_mutation_id in (%s)" % batchstring
			ret = search_db(cursor,qry, verbose=False)

########################################
def decorate_mutations(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:

		time0 = time.time()
		print("====================")
		print("using annotations from  ", table, "to add reliability info to  mutations; pid:", os.getpid())
		update_reliability_column_in_mutation_tables(cursor, table)
		time1 = time.time()
		print(("\t\t %s done in %.3f mins" % (table, float(time1-time0)/60)), os.getpid())

	cursor.close()
	db.close()

	return

#########################################
#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")

	########################
	# set reliability info to 0 in mutation tables - in case we were mucking around with it already
	print("setting reliability to 0 in all mutatitions_chrom_% tables ...")
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like 'mutations_chrom%'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	for table in tables:
		qry = "update %s set reliability_estimate=0" % table
		search_db(cursor,qry)
	print("... done")

	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	table_size = get_table_size(cursor, 'icgc', tables)
	cursor.close()
	db.close()
	# parallelize on those
	tables_sorted = sorted(tables, key=lambda t: table_size[t], reverse=True)
	shuffle(tables_sorted)
	number_of_chunks = 8
	print("number of pll chunks", number_of_chunks)
	parallelize(number_of_chunks, decorate_mutations, tables_sorted, [], round_robin=True)


	return



#########################################
if __name__ == '__main__':
	main()
