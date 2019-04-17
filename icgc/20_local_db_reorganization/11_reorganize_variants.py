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


#########################################
def insert (cursor, table, columns, values):

	nonempty_values = []
	corresponding_columns = []
	for i in range(len(values)):
		if not values[i] or  values[i] == "": continue
		nonempty_values.append(values[i])
		corresponding_columns.append(columns[i])
	qry = "insert into %s (%s) " %(table, ",".join(corresponding_columns))
	qry += "values (%s) " % ",".join(nonempty_values)
	search_db(cursor, qry)


#########################################
def reorganize_donor_variants(cursor, table, columns):

	variants_table = table.replace("_temp","")
	donors = get_donors(cursor, table)
	for donor in donors:
		variants  = set([])
		total_entries = 0
		qry  = "select * from %s where icgc_donor_id='%s' " % (table, donor)
		qry += "and gene_affected is not null and gene_affected !='' "
		ret  = search_db (cursor, qry)
		if not ret: continue # it happens, check "select * from ALL_simple_somatic_temp where icgc_donor_id='DO282'"
		for fields in ret:
			total_entries += 1
			named_field = dict(list(zip(columns,fields)))
			variant_values = []
			for name in variant_columns:
				variant_values.append(quotify(named_field[name]))
			variants.add(",".join(variant_values)) # set => getting rid of duplicates

		for variant in variants:
			insert(cursor, variants_table, variant_columns, variant.split(","))


#########################################
def reorganize(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:

		# the tables should all have the same columns
		qry = "select column_name from information_schema.columns where table_name='%s'"%table
		columns = [field[0] for field in  search_db(cursor,qry)]
		# line by line: move id info into new table
		# for mutation and location check if the info exists; if not make new entry
		time0 = time.time()
		print("====================")
		print("reorganizing variants from ", table, os.getpid())
		reorganize_donor_variants(cursor, table, columns)
		time1 = time.time()
		print(("\t\t %s done in %.3f mins" % (table, float(time1-time0)/60)), os.getpid())

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

	number_of_chunks = 14  # myISAM does not deadlock
	tables_sorted = sorted(tables, key=lambda t: table_size[t], reverse=True)
	parallelize(number_of_chunks, reorganize, tables_sorted, [], round_robin=True)



#########################################
if __name__ == '__main__':
	main()
