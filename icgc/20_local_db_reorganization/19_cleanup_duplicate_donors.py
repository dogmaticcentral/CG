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
from config import Config
from icgc_utils.common_queries  import  *

def count_entries(cursor, somatic_table, icgc_specimen_id):
	qry = "select count(*) from {} where icgc_specimen_id='{}' ".format(somatic_table, icgc_specimen_id)
	ret = search_db(cursor,qry)
	if not ret or ret[0][0]==0: return 0
	return ret[0][0]

#########################################
#########################################
def main():

	#print("disabled - this script deletes certain rows ") # comment out to run
	#exit(1)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which  somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_donor'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,"icgc")
	for table in tables:
		tumor = table.split("_")[0]
		somatic_table =  "%s_simple_somatic" % tumor
		qry  = "select submitted_donor_id, count(*) as c "
		qry += "from %s  group by submitted_donor_id having c>1 " % table
		ret  = search_db(cursor,qry)
		if not ret:
			#print("no duplicates found")
			continue
		print("\n====================")
		print("%s has %d duplicates" % (table,len(ret)))
		for line in ret:
			print("\t", line)
			[submitted_donor_id, count] = line
			qry = "select icgc_donor_id from %s where submitted_donor_id = '%s' " % (table, submitted_donor_id)
			icgc_donor_ids = [ret2[0] for ret2 in search_db(cursor,qry)]
			qry  = "select icgc_specimen_id, icgc_donor_id, specimen_type, tumour_histological_type "
			qry += "from %s_specimen where icgc_donor_id in (%s)" %(tumor, ",".join(["'%s'"%id for id in icgc_donor_ids]))
			ret3 = search_db(cursor,qry)
			specimens = {}
			for line in ret3:
				[icgc_specimen_id, icgc_donor_id, specimen_type, tumour_histological_type] = line
				if not specimen_type in specimens: specimens[specimen_type] = []
				specimens[specimen_type].append([icgc_specimen_id, icgc_donor_id])
			max_mutations = -1
			max_mutation_spec_id = ""
			all_spec_ids = set()
			for spectype, ids in specimens.items():
				for icgc_specimen_id, icgc_donor_id in ids:
					number_of_entries = count_entries(cursor, somatic_table, icgc_specimen_id)
					if number_of_entries==0: continue
					print("\t\t", spectype, icgc_specimen_id, icgc_donor_id, number_of_entries)
					if spectype=='Primary tumour - solid tissue':
						all_spec_ids.add(icgc_specimen_id)
						if max_mutations<number_of_entries:
							max_mutations = number_of_entries
							max_mutation_spec_id = icgc_specimen_id
					else:
						qry = "delete from {} where icgc_specimen_id='{}'".format(somatic_table, icgc_specimen_id)
						search_db(cursor,qry)
			# this should refer to 'Primary tumour - solid tissue' only
			if len(all_spec_ids)<2: continue
			all_spec_ids.remove(max_mutation_spec_id)
			print("\t\t", max_mutations, "for", max_mutation_spec_id, " removing", all_spec_ids)
			for icgc_specimen_id in all_spec_ids:
				qry = "delete from {} where icgc_specimen_id='{}'".format(somatic_table, icgc_specimen_id)
				search_db(cursor,qry)
	cursor.close()
	db.close()


	return



#########################################
if __name__ == '__main__':
	main()
