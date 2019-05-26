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

def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, "icgc")

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]

	for table in tables:
		print(table)
		tumor_short = table.split("_")[0]
		# mut_count = mutation_count_per_donor(cursor, table)
		mut_count = genes_per_patient_breakdown(cursor, table)
		for donor, ct in mut_count.items():
			if ct< 1000: continue
			print(donor, ct)
			qry = "update %s set reliability_estimate=0 where icgc_donor_id='%s'" % (table, donor)
			error_intolerant_search(cursor, qry)
		print()


	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
