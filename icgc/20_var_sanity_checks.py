#! /usr/bin/python

from icgc_utils.common_queries   import  *

def main():

	db     = connect_to_mysql()
	cursor = db.cursor()
	switch_to_db(cursor, "icgc")

	approved_symbol  = 'RPL5'
	rpl5_muts  = mutations_in_gene(cursor, approved_symbol)

	for old_id in ['BP-4331-01', 'BP-4967-01', 'B0-5400-01']:
		qry = "select distinct icgc_mutation_id, icgc_donor_id  from KIRC_simple_somatic "
		qry += "where submitted_sample_id like '%%%s%%'" % old_id
		kirc_muts = [ret[0] for ret in search_db(cursor,qry)]
		new_ids = [ret[1] for ret in search_db(cursor,qry)]

		rpl5_muts_in_kirc = list(set(kirc_muts)&set(rpl5_muts))
		for mut_id in rpl5_muts_in_kirc:
			print old_id, mut_id, new_ids[rpl5_muts_in_kirc.index(mut_id)]
		'''
			this should return
			BP-4331-01 MU617096 DO19982
			BP-4967-01 MU622082 DO16795
			B0-5400-01 MU620539 DO19308		
		'''



	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
