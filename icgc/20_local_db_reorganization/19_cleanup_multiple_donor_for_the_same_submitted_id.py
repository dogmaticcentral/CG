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

# Resolving: multiple icgc_donor_ids for the same submitter id


from config import Config
from icgc_utils.common_queries  import  *

def count_entries(cursor, somatic_table, icgc_specimen_id):
	qry = "select count(*) from {} where icgc_specimen_id='{}' ".format(somatic_table, icgc_specimen_id)
	ret = search_db(cursor,qry)
	if not ret or ret[0][0]==0: return 0
	return ret[0][0]

#########################################
def duplicates_in_donor_table(cursor, table):
	donor_ids_with_duplicates = []
	qry  = "select submitted_donor_id, count(*) as c "
	qry += "from %s  group by submitted_donor_id having c>1 " % table
	ret  = search_db(cursor,qry)
	if ret: donor_ids_with_duplicates = [r[0] for r in ret]
	return donor_ids_with_duplicates

#########################################
def duplicates_in_variants_table(cursor, donor_table, variants_table, donor_ids_with_duplicates_in_donor_table):
	# variants table does not have submitted_donor_id, only submitted sample id
	duplicates = {}
	for submitted_donor_id in donor_ids_with_duplicates_in_donor_table:
		qry =  "select distinct(icgc_donor_id) from %s  " % donor_table
		qry += "where submitted_donor_id='%s' " % submitted_donor_id
		ret  = search_db(cursor,qry, verbose=False)
		if not ret or len(ret)<=1: continue
		duplicates[submitted_donor_id] = [r[0] for r in ret]
	duplicates_in_variants_table = {}
	for submitted_donor_id, icgc_donor_ids in duplicates.items():
		for icgc_donor_id in icgc_donor_ids:
			qry = "select * from %s  " % variants_table
			qry += "where icgc_donor_id = '%s'  limit 1" % icgc_donor_id
			ret  = search_db(cursor,qry, verbose=False)
			if not ret: continue
			if not submitted_donor_id in duplicates_in_variants_table:
				duplicates_in_variants_table[submitted_donor_id] = []
			duplicates_in_variants_table[submitted_donor_id].append(icgc_donor_id)

	for submitted_donor_id  in duplicates.keys():
		dups = []
		if submitted_donor_id in duplicates_in_variants_table:
			dups = duplicates_in_variants_table[submitted_donor_id]

		if len(dups)<=1:
			duplicates_in_variants_table.pop(submitted_donor_id)
			#print("%s has duplicates in donor table but no duplicates in  the variants table:"%submitted_donor_id)
		else:
			#print("%s has duplicates in the variants table:" % submitted_donor_id, dups)
			pass

	#print (len(duplicates), len(duplicates_in_variants_table))

	return duplicates_in_variants_table


#########################################
#########################################
def main():

	#print("disabled - this script deletes certain rows ") # comment out to run
	#exit(1)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which  donor tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_donor'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,"icgc")

	for donor_table in tables:
		print("\n====================")
		print(donor_table)
		tumor = donor_table.split("_")[0]
		variants_table =  "%s_simple_somatic" % tumor
		specimen_table = "%s_specimen" % tumor

		donor_ids_with_duplicates_in_donor_table = duplicates_in_donor_table(cursor, donor_table)
		if len(donor_ids_with_duplicates_in_donor_table)==0:
			#print("\t no duplicates in donor_table", donor_table)
			continue
		# there are duplicates in the donor table;
		# donor table is just a bookkeeping device - do the duplicates really appear in the variants table?
		dups = duplicates_in_variants_table(cursor, donor_table, variants_table, donor_ids_with_duplicates_in_donor_table)
		if len(dups)==0:
			print("\t no duplicates in variants_table", variants_table)
			continue
		print("%s has %d submitted_donor_id ids with duplicate icgc_donor_ids" % (variants_table,len(dups)))

		for submitted_donor_id, icgc_donor_ids in dups.items():
			print("\t submitted_id: %s   icgc_donor_ids: "%submitted_donor_id, icgc_donor_ids)

			qry  = "select icgc_specimen_id, specimen_type from %s " % specimen_table
			qry += "where icgc_donor_id in (%s)" % ",".join([quotify(id) for id in icgc_donor_ids])
			ret3 = search_db(cursor,qry)
			primary_tumors = set([])
			other_tumors = set([])
			for line in ret3:
				[icgc_specimen_id, spectype] = line
				if 'primary' in spectype.lower():
					primary_tumors.add(icgc_specimen_id)
				else:
					other_tumors.add(icgc_specimen_id)

			print("primary:", primary_tumors)
			print("other:", other_tumors)
			removable_specimen_ids = set([])
			if len(primary_tumors)>0:
				usable_specimen_ids = primary_tumors
				removable_specimen_ids = other_tumors
			elif len(other_tumors)>0:
				usable_specimen_ids = other_tumors
			else:
				continue
			print("usable:", usable_specimen_ids)
			print("removable:", removable_specimen_ids)

			max_mutations = -1
			max_mutation_spec_id = ""
			for icgc_specimen_id in usable_specimen_ids:
				number_of_entries = count_entries(cursor, variants_table, icgc_specimen_id)
				print("\t", icgc_specimen_id, number_of_entries)
				if number_of_entries==0: continue
				if max_mutations<number_of_entries:
					max_mutations = number_of_entries
					max_mutation_spec_id = icgc_specimen_id
			if max_mutations<=0: continue
			usable_specimen_ids.remove(max_mutation_spec_id)
			removable_specimen_ids |= usable_specimen_ids
			print(max_mutation_spec_id)
			print(removable_specimen_ids)
			removable_string = ",".join([quotify(s) for s in removable_specimen_ids])
			qry = "delete from {} where icgc_specimen_id in ({})".format(variants_table, removable_string)
			print(qry)
			print()

			search_db(cursor,qry)

	cursor.close()
	db.close()


	return



#########################################
if __name__ == '__main__':
	main()
