#!/usr/bin/python

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


import os.path
import re, commands
from tcga_utils.mysql import  *
from tcga_utils.utils import  get_expected_fields, is_useful, make_named_fields
from time import time


#########################################
def main():

	print "don't do this - there is a lot of screwed up annotation in some cancers"
	print "while the underlying sequencing might be OK"
	print "(it looks like some of the annotation at least might be referring to different splicing)"
	exit()

	db     = connect_to_mysql()
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='tcga' and table_name like '%_somatic_mutations'"
	tables = [field[0] for field in search_db(cursor,qry)]

	switch_to_db (cursor, "tcga")

	conflicts = {}
	total = {}
	for table in tables:
		qry = "select count(1) from %s " % table
		rows = search_db(cursor, qry)
		if not rows or rows[0][0] ==0: continue
		total[table] = rows[0][0]
		qry = "select count(1) from %s where conflict is not null" % table
		rows = search_db(cursor, qry)
		conflicts[table] = int(rows[0][0])


	tables_sorted = sorted(conflicts, key=conflicts.__getitem__, reverse=True)
	for table in tables_sorted:
		print table, conflicts[table], "conflicts, out of", total[table]
		continue

		expected_fields = get_expected_fields(cursor, table, table)
		# get all conflicted groups
		qry = "select * from %s where conflict is not null" % table
		rows = search_db(cursor, qry)
		if not rows: continue

		existing_fields_by_database_id = dict(zip(map(lambda x: int(x[0]), rows),
												  map(lambda x: make_named_fields(expected_fields, x[1:]), rows)))

		bags = []
		keyset = set (existing_fields_by_database_id.keys())
		for db_id, fields in existing_fields_by_database_id.iteritems():
			conflict_ids = set([int(a.split(" with ")[-1]) for a in fields['conflict'].split(";")])
			if not conflict_ids <= keyset:
				print "error - conflicting ids not reciprocally labelled or not in the database"
				exit(1)
			new_bag = set(conflict_ids) | set([db_id])  # set unioin
			new_bags = []
			for bag in bags:
				if not bag & new_bag: # intersection is empty
					new_bags.append(bag)
				else:
					new_bag |= bag  # add the elements of this bag to the new bag
			new_bags.append(new_bag)
			bags = new_bags

		#ef = existing_fields_by_database_id # I neeed a shorthand
		for bag in bags: # bag is a collection of conflicting ids
			conflicting_ids = list(bag)
			print conflicting_ids[0], conflicting_ids[1:]
			keep_id = int(conflicting_ids[0])
			qry = "update %s set aa_change=NULL, conflict=null where id=%d" % (table, keep_id)
			search_db(cursor, qry)
			for delete_id in conflicting_ids[1:]:
				qry = "delete from %s where  id=%d" % (table, int(delete_id))
				search_db(cursor, qry)


	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
