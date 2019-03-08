#! /usr/bin/python3
from config import Config
from icgc_utils.common_queries  import  *



#########################################
#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which  somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_donor'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,"icgc")
	for table in tables:
		print("\n====================")
		print("inspecting", table)
		tumor = table.split("_")[0]
		somatic_table =  "%s_simple_somatic" % tumor
		qry  = "select submitted_donor_id, count(*) as c "
		qry += "from %s  group by submitted_donor_id having c>1 " % table
		ret  = search_db(cursor,qry)
		if not ret:
			print("no duplicates found")
			continue
		print("has %d duplicates" % len(ret))
		for line in ret:
			print("\t", line)
			[submitted_donor_id, count] = line
			qry = "select icgc_donor_id from %s where submitted_donor_id = '%s' " % (table, submitted_donor_id)
			icgc_donor_ids = [ret2[0] for ret2 in search_db(cursor,qry)]
			qry  = "select icgc_specimen_id, icgc_donor_id, specimen_type, tumour_histological_type "
			qry += "from %s_specimen where icgc_donor_id in (%s)" %(tumor, ",".join(["'%s'"%id for id in icgc_donor_ids]))
			ret3 = search_db(cursor,qry)
			for line in ret3:
				print("\t\t", line)


	cursor.close()
	db.close()


	return



#########################################
if __name__ == '__main__':
	main()
