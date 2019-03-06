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
		print("has duplicates")
		for line in ret:
			print("\t", line)
			[submitted_donor_id, count] = line
			qry = "select icgc_donor_id from %s where submitted_donor_id = '%s' " % (table, submitted_donor_id)
			icgc_donor_ids = [ret2[0] for ret2 in search_db(cursor,qry)]
			qry = "select t1.icgc_mutation_id, t1.icgc_specimen_id, t2.icgc_specimen_id, t1.icgc_sample_id, t2.icgc_sample_id "
			qry += "from %s as t1, %s as t2 " % (somatic_table, somatic_table)
			qry += "where  t1.icgc_donor_id='%s' and t2.icgc_donor_id='%s' " % tuple(icgc_donor_ids)
			qry += "and t1.icgc_mutation_id = t2.icgc_mutation_id"
			print(qry)
			ret3 = search_db(cursor,qry)
			print(ret3[0])
			exit()


	cursor.close()
	db.close()


	return



#########################################
if __name__ == '__main__':
	main()
