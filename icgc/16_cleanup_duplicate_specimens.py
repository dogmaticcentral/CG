#! /usr/bin/python3
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
	qry += "where table_schema='icgc' and table_name like '%_specimen'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,"icgc")
	for table in tables:
		qry  = "select icgc_donor_id, count(*) as c "
		qry += "from %s  group by icgc_donor_id having c>1 " % table
		ret  = search_db(cursor,qry)
		if not ret: # looks like we alwas have duplicates
			print("%s no duplicates found"% table)
			continue
		print("\n====================")
		print("%s has %d duplicates" % (table,len(ret)))

		# do we have duplicates in the variants table?
		tumor = table.split("_")[0]
		somatic_table =  "%s_simple_somatic" % tumor
		icgc_donor_ids = [r[0] for r in ret]
		really_problematic = []
		for icgc_donor_id in icgc_donor_ids:
			qry  = "select distinct(icgc_specimen_id) "
			qry += "from %s where icgc_donor_id='%s' " % (somatic_table,icgc_donor_id)
			ret  = search_db(cursor,qry)
			if ret and len(ret)>1:
				really_problematic.append(icgc_donor_id)
		print(really_problematic)
		print("duplicates in the variants table:", len(really_problematic))
		exit()
		continue
		# most of these are inocuous, with normal sample not appearing in the variants table
		# this, however is not always the case
		for line in ret:
			#print("\t", line)
			[icgc_donor_id, count] = line
			qry = "select icgc_specimen_id from %s where icgc_donor_id = '%s' " % (table, icgc_donor_id)
			icgc_specimen_ids = [ret2[0] for ret2 in search_db(cursor,qry)]
			qry  = "select icgc_specimen_id, icgc_donor_id, specimen_type, tumour_histological_type "
			qry += "from %s_specimen where icgc_specimen_id in (%s)" %(tumor, ",".join(["'%s'"%id for id in icgc_specimen_ids]))
			ret3 = search_db(cursor,qry)
			other_spec_ids = set()
			primary_spec_ids = set()
			normal_spec_ids  = set()
			for line in ret3:
				[icgc_specimen_id, icgc_donor_id, specimen_type, tumour_histological_type] = line
				if 'Primary' in specimen_type:
					primary_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))
				elif 'Normal' in specimen_type:
					normal_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))
				else:
					other_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))
			if len(primary_spec_ids)==0:
				print ("\t no primary ids", icgc_donor_id)
				print ("\t other spec ids", other_spec_ids)
			elif len(primary_spec_ids)>1:
				print ("\t multiple primary ids", icgc_donor_id)
				print ("\t primary spec ids", primary_spec_ids)

		# 		if not specimen_type in specimens: specimens[specimen_type] = []
		# 		specimens[specimen_type].append([icgc_specimen_id, icgc_donor_id])
		# 	max_mutations = -1
		# 	max_mutation_spec_id = ""
		# 	all_spec_ids = set()
		# 	for spectype, ids in specimens.items():
		# 		for icgc_specimen_id, icgc_donor_id in ids:
		# 			number_of_entries = count_entries(cursor, somatic_table, icgc_specimen_id)
		# 			if number_of_entries==0: continue
		# 			print("\t\t", spectype, icgc_specimen_id, icgc_donor_id, number_of_entries)
		# 			if spectype=='Primary tumour - solid tissue':
		# 				all_spec_ids.add(icgc_specimen_id)
		# 				if max_mutations<number_of_entries:
		# 					max_mutations = number_of_entries
		# 					max_mutation_spec_id = icgc_specimen_id
		# 			else:
		# 				qry = "delete from {} where icgc_specimen_id='{}'".format(somatic_table, icgc_specimen_id)
		# 				search_db(cursor,qry)
		# 	# this should refer to 'Primary tumour - solid tissue' only
		# 	if len(all_spec_ids)<2: continue
		# 	all_spec_ids.remove(max_mutation_spec_id)
		# 	print("\t\t", max_mutations, "for", max_mutation_spec_id, " removing", all_spec_ids)
		# 	for icgc_specimen_id in all_spec_ids:
		# 		qry = "delete from {} where icgc_specimen_id='{}'".format(somatic_table, icgc_specimen_id)
		# 		search_db(cursor,qry)
		#exit()
	cursor.close()
	db.close()


	return



#########################################
if __name__ == '__main__':
	main()
