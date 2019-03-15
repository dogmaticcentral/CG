#! /usr/bin/python3
from config import Config
from icgc_utils.common_queries  import  *

def count_entries(cursor, somatic_table, icgc_specimen_id):
	qry = "select count(*) from {} where icgc_specimen_id='{}' ".format(somatic_table, icgc_specimen_id)
	ret = search_db(cursor,qry)
	if not ret or ret[0][0]==0: return 0
	return ret[0][0]

def find_spec_id_with_max_entries(spec_ids_w_description, entries_per_specimen):
	max_count = 0
	max_spec_id = None
	for psi in spec_ids_w_description:
		[spec_id, description] = psi.split(":")
		print ("\t\t ", spec_id, entries_per_specimen[spec_id], description)
		if max_count<entries_per_specimen[spec_id]:
			max_count=entries_per_specimen[spec_id]
			max_spec_id = spec_id
	return max_spec_id

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
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,"icgc")
	for somatic_table in tables:
		tumor = somatic_table.split("_")[0]
		icgc_donor_ids = [r[0] for r in search_db(cursor, "select distinct icgc_donor_id from %s"%somatic_table)]
		problematic = []
		specimen_ids = {}
		for icgc_donor_id in icgc_donor_ids:
			qry  = "select distinct(icgc_specimen_id) "
			qry += "from %s where icgc_donor_id='%s' " % (somatic_table,icgc_donor_id)
			ret  = search_db(cursor,qry)
			if ret and len(ret)>1:
				problematic.append(icgc_donor_id)
				specimen_ids[icgc_donor_id] = [r[0] for r in ret]
		if len(problematic)==0:
			print("%s has no duplicates in the variant table"% somatic_table)
			continue

		print("%s has %d donor ids with duplicate specimen ids in the variant table" % (somatic_table,len(problematic)))

		# most of these are innocuous, with normal sample not appearing in the variants table
		# this, however is not always the case
		for icgc_donor_id in problematic:
			print ("specimen ids for %s" % icgc_donor_id, specimen_ids[icgc_donor_id])

			qry  = "select icgc_specimen_id, count(*) as c  from %s " % somatic_table
			qry += "where icgc_donor_id = '%s' and reliability_estimate=1 " % icgc_donor_id
			qry += "group by  icgc_specimen_id"
			ret2 = search_db(cursor,qry)
			entries_per_specimen = dict(ret2)
			print(entries_per_specimen)
			qry  = "select icgc_specimen_id, icgc_donor_id, specimen_type, tumour_histological_type "
			qry += "from %s_specimen where icgc_specimen_id in (%s)" % (tumor, ",".join(["'%s'"%id for id in entries_per_specimen.keys()]))
			ret3 = search_db(cursor,qry)
			if not ret3:
				search_db(cursor,qry,verbose=True)
				exit(1)
			other_spec_ids = set()
			primary_spec_ids = set()
			normal_spec_ids  = set()
			metastatic_spec_ids  = set()
			removable_ids = set(specimen_ids[icgc_donor_id])
			for line in ret3:
				[icgc_specimen_id, icgc_donor_id, specimen_type, tumour_histological_type] = line
				if 'Primary' in specimen_type:
					primary_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))
				elif 'Normal' in specimen_type:
					normal_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))
				elif 'Metastatic' in specimen_type:
					metastatic_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))
				else:
					other_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))

			#################################
			keep_id = None
			if len(primary_spec_ids)==0:
				print ("\t no primary ids", icgc_donor_id)
				print ("\t other spec ids", other_spec_ids)
				if len(metastatic_spec_ids)==1:
					metastatic_spec_id = metastatic_spec_ids.pop().split(":")[0]
					print ("there is only one metastastic tumor id: ", metastatic_spec_id)
					keep_id = metastatic_spec_id

				elif len(metastatic_spec_ids)>1:
					max_spec_id = find_spec_id_with_max_entries(metastatic_spec_ids, entries_per_specimen)
					if max_spec_id:
						keep_id =  max_spec_id
						print ("using  metastastic tumor id: ", max_spec_id)
				else: # if all else fails, use whatever specimens are available
					max_spec_id = find_spec_id_with_max_entries(other_spec_ids, entries_per_specimen)
					if max_spec_id:
						keep_id =  max_spec_id
						print ("using  'other' tumor id: ", max_spec_id)

			elif len(primary_spec_ids)>1:
				print ("\t multiple primary ids", icgc_donor_id)
				max_spec_id = find_spec_id_with_max_entries(primary_spec_ids, entries_per_specimen)
				if max_spec_id:
					keep_id =  max_spec_id
			else:
				primary_spec_id = primary_spec_ids.pop().split(":")[0]
				print ("there is only one primary id: ", primary_spec_id)
				keep_id = primary_spec_id

			##
			if keep_id:
				removable_ids.remove(keep_id)
				removable_ids_string = ",".join(["'%s'"%rid for rid in removable_ids])
				qry = "delete from {} where icgc_specimen_id in ({})".format(somatic_table, removable_ids_string)
				print(qry)
				search_db(cursor,qry)

	cursor.close()
	db.close()


	return



#########################################
if __name__ == '__main__':
	main()
