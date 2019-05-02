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

variant_columns = ['icgc_mutation_id', 'chromosome','icgc_donor_id', 'icgc_specimen_id', 'icgc_sample_id',
				   'submitted_sample_id','control_genotype', 'tumor_genotype', 'total_read_count', 'mutant_allele_read_count']


def check_and_drop(cursor, table):
	if check_table_exists(cursor, 'icgc', table):
		search_db(cursor, "drop table %s"% table)
	return

def time_qry(cursor, qry):
	time0 = time.time()
	search_db(cursor,qry, verbose=False)
	time1 = time.time()
	print("\n%s\ndone in %.3f mins" % (qry, float(time1-time0)/60))

#########################################
def reorganize_variants(cursor, orig_icgc_table):

	keep_string = ", ".join([c for c in variant_columns])

	switch_to_db(cursor, "icgc")
	tmp_table = "scratch_%d_%s" % (os.getpid(), orig_icgc_table)
	check_and_drop(cursor, tmp_table)
	qry = "create temporary table  %s  engine=MYISAM  as select %s from %s " % (tmp_table, keep_string, orig_icgc_table)
	time_qry(cursor,qry)

	new_icgc_table = orig_icgc_table.replace("_temp", "")
	check_and_drop(cursor,new_icgc_table)
	qry = "create table %s like %s" % (new_icgc_table, tmp_table)
	time_qry(cursor,qry)

	qry = "insert into %s select distinct * from %s" % (new_icgc_table, tmp_table)
	time_qry(cursor,qry)

	# add back the primary key
	qry = "alter table %s  add column id  int not null primary key auto_increment first" % new_icgc_table
	time_qry(cursor,qry)

	check_and_drop(cursor,tmp_table)


#########################################
def reorganize(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:

		time0 = time.time()
		print("====================")
		print("reorganizing variants from ", table, os.getpid())
		reorganize_variants(cursor, table)
		# TODO: add mut_to_total_read_count_ratiom pathogenicity and reliability columns to the new variants table
		time1 = time.time()
		print(("\t\t %s (%d) done in %.3f mins" % (table, tables.index(table),  float(time1-time0)/60)), os.getpid())

	cursor.close()
	db.close()

	return

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
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	table_size = get_table_size(cursor, 'icgc', tables)
	cursor.close()
	db.close()

	tables_sorted = sorted(tables, key=lambda t: table_size[t], reverse=True)
	half = int(len(tables_sorted)/2)
	tables_mirrored  = tables_sorted[0:half] + list(reversed(tables_sorted[half:]))
	number_of_chunks = half

	parallelize(number_of_chunks, reorganize, tables_mirrored, [], round_robin=True)



#########################################
if __name__ == '__main__':
	main()
