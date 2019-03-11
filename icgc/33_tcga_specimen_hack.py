#! /usr/bin/python3

# tcga does not have specim columns
# so we will take the TCGA submitted sample id to play the role
# we can only hope we do not have multiple samples from the same specimen
# bcs we do not know what happens then

from icgc_utils.common_queries   import  *
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")

	for table in tables:
		print(table)
		qry  = "update %s " % table
		qry += "set icgc_specimen_id=submitted_sample_id "
		qry += "where icgc_donor_id like 'DOT_%'"
		print(qry)
		search_db(cursor,qry)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
