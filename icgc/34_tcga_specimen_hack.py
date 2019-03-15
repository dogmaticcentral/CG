#! /usr/bin/python3

# tcga does not have specim columns
# so we will take the TCGA submitted sample id to play the role
# we can only hope we do not have multiple samples from the same specimen
# bcs we do not know what happens then

from icgc_utils.common_queries   import  *
from icgc_utils.icgc   import  *
from config import Config

def spec_from_TCGA(sample_barcode):
	# we want to translate this to something similar to what ICGC is using
	# roughly: Normal, Primary, Metastatic, Recurrent, Cell_line
	tcga_sample_code = sample_barcode.split("-")[3][:2]
	if tcga_sample_code in ['01','03','05', '09']:
		return 'Primary'

	elif tcga_sample_code in ['02','04','40']:
		return 'Recurrent'

	elif tcga_sample_code in ['06','07']:
		return 'Metastatic'

	elif tcga_sample_code in ['10','11','12','13','14']:
		return 'Normal'

	elif tcga_sample_code in ['08','50']:
		return 'Cell_line'

	return "Other"


def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	somatic_tables = [field[0] for field in search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")

	for somatic_table in somatic_tables:
		print(somatic_table)
		qry  = "update %s " % somatic_table
		qry += "set icgc_specimen_id=submitted_sample_id "
		qry += "where icgc_donor_id like 'DOT_%'"
		#print(qry)
		search_db(cursor,qry)
		# do we have the corresponding specimen table?
		specimen_table = somatic_table.replace("simple_somatic", "specimen")
		if not check_table_exists(cursor,'icgc', specimen_table):
			make_specimen_table(cursor, 'icgc', specimen_table)
			print(specimen_table)
		qry  = "select distinct(d.icgc_donor_id), s.icgc_specimen_id  "
		qry += "from %s d inner join %s s " % (somatic_table, somatic_table)
		qry += "on  d.icgc_donor_id=s.icgc_donor_id where  d.icgc_donor_id like 'DOT_%'"
		ret = search_db(cursor,qry)
		if not ret: continue
		if 'Error' in ret[0][0]:
			search_db(cursor,qry,verbose=True)
			exit(1)
		for icgc_donor_id, icgc_specimen_id in ret:
			specimen_type = spec_from_TCGA(icgc_specimen_id)
			#print(icgc_donor_id, icgc_specimen_id, specimen_type)
			# check specimen registered in specimen table
			qry = "select icgc_donor_id, icgc_specimen_id from %s " % specimen_table
			qry += "where icgc_specimen_id='%s' " % icgc_specimen_id
			ret = search_db(cursor,qry)
			if ret:
				#print(ret)
				continue
			# if it does not exist store as new
			qry = "select id from %s order by id desc limit 1 "  %  specimen_table
			ret = search_db(cursor,qry)
			last_id_in_specimen_table = 0 if not ret else ret[0][0]
			new_id = last_id_in_specimen_table + 1
			qry = "insert into %s (id, icgc_specimen_id, icgc_donor_id, specimen_type) " %  specimen_table
			qry += "values ('{}','{}','{}','{}')".format(new_id, icgc_specimen_id, icgc_donor_id, specimen_type)
			search_db(cursor,qry)
			#exit()
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
