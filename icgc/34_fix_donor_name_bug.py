#! /usr/bin/python

# tcga does not have specimen id columns
# so we will take the TCGA submitted sample id to play the role
# we can only hope we do not have multiple samples from the same specimen
# bcs we do not know waht happens then

from icgc_utils.common_queries   import  *
from icgc_utils.processes import *

def fix(tables, other_args):

	db     = connect_to_mysql()
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")

	for table in tables:
		print(table)
		if True:
			qry  = "select distinct(SUBSTRING(submitted_sample_id,1,12)) "
			qry += "from %s " % table
			qry += "where  submitted_sample_id like 'TCGA%'"
		#else:
		#	qry  = "select distinct(SUBSTRING(submitted_sample_id,1,16)) "
		#	qry += "from %s " % table
		#	qry += "where  submitted_sample_id like 'TARGET%'"
		ret = search_db(cursor,qry)
		if not ret: continue
		tcga_donors = [r[0] for r in ret]
		print('tcga_donors:', len(tcga_donors))
		count = 0
		for td in  tcga_donors:
			if tcga_donors.index(td)%100 == 0: print("\t", tcga_donors.index(td), "donors")
			qry  = "select distinct(icgc_donor_id) "
			qry += "from %s " % table
			qry += "where submitted_sample_id like '%s%%' " %  td
			ret = search_db(cursor,qry)
			if not ret: continue
			icgc_donor_ids = [r[0] for r in ret]
			if len(icgc_donor_ids)<2: continue
			print(td, icgc_donor_ids)
			count += 1
			if count%10==0: print(table, count)
			dos  = [d for d in icgc_donor_ids if d[2]!='T']
			dots = [d for d in icgc_donor_ids if d[2]=='T']
			if len (dos) == 0:
				print(table, td, "no original ids")
				dos = dots[:1]
				dots = dots[1:]
			if len(dos)>1:
				print(table, td, "too many original ids")
				print(dos)
				continue
			do = dos[0]

			for dot in dots:
				mut_ids = [r for r in  search_db(cursor, "select id, icgc_mutation_id from %s where icgc_donor_id='%s'" % (table,dot))]
				for new_variant_id, mut_id in mut_ids:
					if mut_id[2]=='T':
						# this came from tcga -- just change the donor id
						change_donor_id = True
						pass
					else:
						qry  = "select id from %s where icgc_donor_id='%s' " % (table,do)
						qry += "and icgc_mutation_id='%s'" % mut_id
						ret  = search_db(cursor,qry)
						if ret: # the old variant entry exists
							change_donor_id = False
						else:
							change_donor_id = True

					if change_donor_id:
						qry = "update %s set icgc_donor_id = '%s' where id=%d " % (table, do, new_variant_id)
					else: # drop entry
						qry = "delete from %s where id=%d " % (table, new_variant_id)

					search_db(cursor,qry)

		print(count)
	cursor.close()
	db.close()


#########################################
#########################################
def main():



	db     = connect_to_mysql()
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	number_of_chunks = 10 # 8  # myISAM does not deadlock
	parallelize(number_of_chunks, fix, tables, [])




#########################################
if __name__ == '__main__':
	main()
