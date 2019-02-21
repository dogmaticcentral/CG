#! /usr/bin/python3

from icgc_utils.mysql import *

#########################################
def main():

	db      = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor  = db.cursor()
	db_name = "homo_sapiens_core_94_38"
	switch_to_db(cursor, db_name)

	qry  = "select t.stable_id, g.stable_id, g.canonical_transcript_id from transcript as t "
	qry += "left join gene as g on t.gene_id=g.gene_id "
	ret = search_db(cursor, qry)
	if ret:
		for [transcript_stable_id, gene_stable_id,  canonical_transcript_id] in ret:
			ret2 =  search_db(cursor, "select stable_id from transcript where transcript_id=%s"%canonical_transcript_id)
			if not ret2: continue
			canonical_stable = ret2[0][0]
			print("\t".join([transcript_stable_id, gene_stable_id, canonical_stable]))
	else:
		print("No ret for", qry)

	cursor.close()
	db.close()

	return True


#########################################
if __name__ == '__main__':
	main()
