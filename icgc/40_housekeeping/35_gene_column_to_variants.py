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


#########################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	somatic_tables = [field[0] for field in hard_landing_search(cursor,qry)]
	for somatic_table in somatic_tables:
		print ("adding gene_symbol column to", somatic_table)
		# this won't do anything if the column exists
		add_column(cursor, 'icgc', somatic_table, "gene_symbol", "VARCHAR(30)", default=None, after_col="chromosome")
		print ("filling in gene_symbol in", somatic_table)
		qry = "update %s as s, mutation2gene as m set s.gene_symbol=m.gene_symbol " % somatic_table
		qry += "where s.icgc_mutation_id = m.icgc_mutation_id"
		error_intolerant_search(cursor, qry)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
