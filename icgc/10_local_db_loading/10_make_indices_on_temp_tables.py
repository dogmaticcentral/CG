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

import time

from config import Config
from icgc_utils.mysql   import  *
from icgc_utils.processes import *

def make_indices(tables, other_args):
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")
	for table in tables:
		print(table)
		time0 = time.time()
		qry  = "create index somatic_mut_idx on %s (icgc_mutation_id)" % table
		search_db(cursor,qry,verbose=True)
		qry  = "create index somatic_donor_idx on %s (icgc_donor_id)" % table
		search_db(cursor,qry,verbose=True)
		qry  = "create index sample_idx on %s (submitted_sample_id)" % table
		search_db(cursor,qry,verbose=True)
		qry  = "create index mut_gene_idx on %s (icgc_mutation_id, gene_affected)" %  table
		search_db(cursor,qry,verbose=True)
		qry  =  "create index chrom_start_pos_idx on %s (chromosome, start_position)" % table
		search_db(cursor, qry, verbose=True)
		qry  =  "create index chrom_end_pos_idx on %s (chromosome, end_position)" % table
		search_db(cursor, qry, verbose=True)
		print(("\t\t  %s done in %.3f mins" % (table, float(time.time()-time0)/60)))
	cursor.close()
	db.close()


#########################################
#########################################
def main():


	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")

	# indices on simple somatic temp
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	number_of_chunks = 1
	parallelize(number_of_chunks, make_indices, tables, [])

#########################################
if __name__ == '__main__':
	main()

