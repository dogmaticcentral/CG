#! /usr/bin/python


from icgc_utils.mysql   import  *

#########################################
#########################################
def main():


	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################

	switch_to_db(cursor,"icgc")
	for table in tables:
		print "================================="
		print table
		# total number of donors?
		qry  = "select distinct(icgc_donor_id) from %s " % table
		donors = [ret[0] for ret in search_db(cursor,qry)]
		print "\t donors: ", len(donors)
		qry  = "select distinct(icgc_specimen_id) from %s " % table
		specimens = [ret[0] for ret in search_db(cursor,qry)]
		print "\t specimens: ", len(specimens)
		#==================
		qry = "select icgc_donor_id, count(distinct(icgc_specimen_id)) as c "
		qry +="from %s group by icgc_donor_id having c>1 " % table
		ret = search_db(cursor,qry)
		if not ret: continue
		for line in ret:
			[icgc_donor_id,count] = line
			qry = "select  distinct(icgc_specimen_id)  from %s " % table
			qry += "where icgc_donor_id = '%s'" % icgc_donor_id
			spec_ids = [r[0] for r in search_db(cursor,qry)]
			#print "\t", icgc_donor_id
			#for spec_id in spec_ids:
			#	qry = " select specimen_type from %s " % table.replace("simple_somatic","specimen")
			#	qry += "where icgc_specimen_id = '%s'" % spec_id
			#	spec_type = search_db(cursor,qry)[0][0]
			#   	print "\t\t", spec_id, spec_type
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
