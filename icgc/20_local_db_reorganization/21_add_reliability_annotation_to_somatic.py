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

from icgc_utils.common_queries import *
from icgc_utils.processes import *
from config import Config


#########################################
def add_columns(cursor, table):

	# mut_to_total_read_count_ratio
	column_name = "reliability_estimate"
	if not column_exists(cursor, "icgc", table, column_name):
		qry = "ALTER TABLE %s ADD  %s boolean  default 0" % (table, column_name)
		search_db(cursor,qry, verbose=True)


#########################################
def add_reliability_annotation(tables, other_args):
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:
		print(table, "running")
		qry  = "update %s set mut_to_total_read_count_ratio=mutant_allele_read_count/total_read_count " % table
		# for some data sets we do not have the read count info
		qry += "where total_read_count is not null and total_read_count>0 "
		qry += "and mutant_allele_read_count is not null"
		search_db(cursor,qry, verbose=False)
		qry  = "update %s set reliability_estimate = 1 " % table
		# for some data sets we do not have the read count info
		qry += "where total_read_count is null or "
		qry += "(mutant_allele_read_count>=10 and mut_to_total_read_count_ratio>=0.2)"
		search_db(cursor,qry, verbose=False)

	cursor.close()
	db.close()


#########################################
#########################################
def main():

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	for table in tables:
		print("checking/adding reliability column to", table)
		add_columns(cursor, table)
	cursor.close()
	db.close()

	print("adding annotation")

	number_of_chunks = 20
	parallelize(number_of_chunks, add_reliability_annotation, tables, [])


	return



#########################################
if __name__ == '__main__':
	main()
