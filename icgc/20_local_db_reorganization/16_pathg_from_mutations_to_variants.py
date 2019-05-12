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
# an attempt to solve this by a mysql one-liner
# didn't work for me, to put it mildly
def fix_pathogenicity(tables, other_args):
	chromosomes = other_args[0]
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table  in tables:
		print()
		print("====================")
		print("processing  ", table, os.getpid())
		for chromosome in chromosomes:
			time0 = time.time()
			print("\t {} chrom {}".format(table, chromosome))
			qry  = "select icgc_mutation_id, id from %s " % table
			qry += "where chromosome='%s'" % chromosome
			ret = search_db(cursor,qry, verbose=False)
			if not ret: continue

			table_ids = {}
			for icgc_mutation_id, id in ret:
				if not icgc_mutation_id in table_ids: table_ids[icgc_mutation_id] = []
				table_ids[icgc_mutation_id].append(id)

			qry  = "select icgc_mutation_id from mutations_chrom_%s " % chromosome
			qry += "where pathogenicity_estimate=1"
			ret = search_db(cursor,qry, verbose=False)
			pathogenic_muts = set([r[0] for r in ret])

			path_muts_in_table = set(table_ids.keys()).intersection(pathogenic_muts)
			path_ids = []
			for p in path_muts_in_table: path_ids.extend(table_ids[p])

			for batch in range(0,len(path_ids),100):
				qry  = "update %s " % table
				qry += "set pathogenicity_estimate=1 where  id in (%s)" % ",".join([str(i) for i in path_ids[batch:batch+100]])
				search_db(cursor,qry, verbose=False)
			print("\t %s, chrom %s done  in %.3f mins" % (table, chromosome, float(time.time()-time0)/60))

	cursor.close()
	db.close()
	return


#########################################
def main():

	chromosomes = [str(i) for i in range(1,22)] + ["X","Y"]

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	table_size = get_table_size(cursor, 'icgc', tables)
	cursor.close()
	db.close()

	tables_sorted = sorted(tables, key=lambda t: table_size[t], reverse=False)
	shuffle(tables_sorted) # so called 'lazy load balancing'
	number_of_chunks = 8

	parallelize(number_of_chunks, fix_pathogenicity,  tables_sorted, [chromosomes], round_robin=True)



#########################################
if __name__ == '__main__':
	main()
