#!/usr/bin/python -u
# store the meta info about the maf files: name and the reference genome,  for now
#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2016 Ivana Mihalek.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#
# store the meta info about the maf files: name and the reference genome,  for now

from old_tcga_tools.tcga_utils.utils import *
from random import random
from old_tcga_tools.tcga_utils.ucsc import segment_from_das

##################################################################################
##################################################################################
def main():

	db = connect_to_mysql()
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='tcga' and table_name like '%_somatic_mutations'"
	tables = [field[0] for field in search_db(cursor,qry)]

	switch_to_db(cursor, "tcga")

	for table in tables:
		# total number of samples:
		qry = "select count(distinct tumor_sample_barcode) from %s" % table
		total_samples  = search_db(cursor,qry)[0][0]
		mutations_meta = table.split("_")[0]+"_mutations_meta"
		qry  = "select * from %s where diagnostics like '%%stutter%%'" % mutations_meta
		rows = search_db(cursor,qry)
		if not rows: continue
		print
		print " ** ", table, "total samples:", total_samples
		for row in rows:
			meta_id = row[0]
			terms = row[-1].split(";")
			for term in terms:
				if not "stutter" in term: continue
				sample_ids = term.split("=>")[-1].replace(" ","").split(",")
				sample_ids_quoted = ",".join(map(lambda x: '"' + x + '"', sample_ids))
				# make sure it is the same meta_info_id,
				# because we might have the replacement
				qry  = "delete from %s where tumor_sample_barcode in (%s) " % (table, sample_ids_quoted)
				qry += "and meta_info_id=%d" % meta_id
				search_db(cursor,qry)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
