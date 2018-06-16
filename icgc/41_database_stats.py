#! /usr/bin/python


from icgc_utils.common_queries   import  *

verbose = False


#########################################
#########################################
# produce table of the format
# tumor short | tumor long | number of patients | avg number of mutations per patient |
#  number of patients with mutated rpl5 (%of patients; number of genes which are seen mutated in the same or bigger number of patients)
#  | ditto for rp111

def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")

	total_donors = 0
	donors_with_multiple_specimens = 0
	for table in tables:
		tumor_short = table.split("_")[0]
		fields = [tumor_short]
		if verbose: print "================================="
		if verbose: print table

		# total number of donors?
		qry    = "select count(distinct icgc_donor_id) from  %s " % table
		donors = search_db(cursor,qry)[0][0]
		total_donors += donors

		# specimens per donor?
		qry  = "select  icgc_donor_id, count(distinct  icgc_sample_id) ct "
		qry += "from  %s  " % table
		qry += "group by icgc_donor_id having ct>1 order by ct desc"
		ret = search_db(cursor,qry)
		if not ret: continue
		donors_with_multiple_specimens += len(ret)

	print "total_donors:", total_donors
	print "donors_with_multiple_specimens:", donors_with_multiple_specimens
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
