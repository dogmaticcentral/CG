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
#
# I've managed to screw this one up (missed splice?), so go back and fix

import subprocess
import time, re

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *
from random import shuffle
from config import Config

# this is set literal
mutation_pathogenic = {'missense','frameshift',  'stop_gained', 'inframe',
			  'stop_lost', 'inframe_deletion', 'inframe_insertion',
			  'start_lost', 'disruptive_inframe_deletion',
			   'exon_loss', 'disruptive_inframe_insertion',
			  'splice', '5_prime_UTR_premature_start_codon_gain',
			  'splice_acceptor', 'splice_region', 'splice_donor'
			 }

location_pathogenic = { 'splice', '5_prime_UTR_premature_start_codon_gain',
			  'splice_acceptor', 'splice_region', 'splice_donor',
}


#########################################
def fix_pathogenicity(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table  in tables:
		print()
		print("====================")
		print("processing  ", table, os.getpid())
		add_boolean_column(cursor, 'icgc', table, 'pathogenicity_estimate')

		qry  = "select id, icgc_mutation_id, chromosome from %s " % table
		for line in search_db(cursor,qry):
			[variant_id, icgc_mutation_id, chromosome] = line
			qry2  = "update %s v, mutations_chrom_%s m " % (table, chromosome)
			qry2 += "set v.pathogenicity_estimate=m.pathogenicity_estimate "
			qry2 += "where v.id = %s " % variant_id
			qry2 += "and  m.icgc_mutation_id = '%s'" % icgc_mutation_id
			search_db(cursor,qry2)

	cursor.close()
	db.close()
	return


#########################################
def main():
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	number_of_chunks = 12  # myISAM does not deadlock
	parallelize(number_of_chunks, fix_pathogenicity, tables, [], round_robin=True)



#########################################
if __name__ == '__main__':
	main()
