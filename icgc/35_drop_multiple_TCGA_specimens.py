#! /usr/bin/python3

# this was the idea:
# I don't know what to do with multiple specimens from the same
# donor - is that the same tumor or not?
# keep them in the following order of preference
# primary, metastatic, recurrent, cell line
# if here are two samples from the same stage,
# choose one from ICGC, if both from tcga choose any

# however,
# it turns out that in many cases the overlap is 0
# so just keep plowing with what I have

from icgc_utils.common_queries  import  *
from config import Config

def cleanup (cursor, table, donor, specimen_ids):
	print("\t", donor)
	mutations = {}
	common_muts = set([])
	for spec_id in specimen_ids:
		spec_type = spec_type_res(cursor, table, spec_id)[0]
		mutations[spec_id]  = [m[0] for m in  search_db(cursor,"select  icgc_mutation_id from %s where icgc_specimen_id='%s'" % (table,spec_id))]
		if len(common_muts)==0:
			common_muts = set(mutations[spec_id] )
		else:
			common_muts &= set(mutations[spec_id] )
		print("\t", spec_id, spec_type, len(mutations[spec_id]), mutations[spec_id][:3])
	print(len(common_muts))


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


def spec_type_res(cursor, table, spec_id):
	tumor_short = table.split("_")[0]
	if spec_id[:2] == 'SP':
		return get_specimen_type(cursor, tumor_short,[spec_id])[spec_id]
	if spec_id[:4]=='TCGA' or spec_id[:4]=='TARG':
		return spec_from_TCGA(spec_id)
	else:
		return "unk"

def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc'  and table_name like '%simple_somatic'"

	tables = [field[0] for field in search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")

	for table in tables:
		print(table)
		donors = get_donors(cursor,table)

		for donor in donors:
			qry  = "select distinct icgc_specimen_id from %s " % table
			qry += "where icgc_donor_id = '%s' " % donor
			specimen_ids = [r[0] for  r in search_db(cursor,qry)]
			if len(specimen_ids)>1:
				print("multiple specimen ids for", donor)
				#cleanup (cursor,table, donor, specimen_ids)


	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
