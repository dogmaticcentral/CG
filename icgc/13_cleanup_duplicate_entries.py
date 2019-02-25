#! /usr/bin/python3


from .icgc_utils.common_queries  import  *



#########################################
#########################################
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,"icgc")
	for table in tables:
		print("====================")
		print("inspecting ", table)
		# coumn names/headers
		colnames = get_column_names (cursor, "icgc", table)
		# base name
		tumor_short  = table.split("_")[0]
		# a hack to get all entries that have all relevant ids identical
		qry = "select concat(icgc_mutation_id,'_', icgc_donor_id,'_',icgc_specimen_id,'_',icgc_sample_id) as mega_id, "
		qry += "count(*) as c from %s  group by mega_id having c>1 " % table
		ret  = search_db(cursor,qry)
		if not ret:
			print("no duplicates found\n")
			continue

		# loop over all duplicate entries
		for line in ret:
			[mega_id, ct] = line
			print(mega_id, ct)
			[icgc_mutation_id, icgc_donor_id, icgc_specimen_id, icgc_sample_id] = mega_id.split("_")
			# donor
			qry  = "select * from %s_donor " % tumor_short
			qry += "where icgc_donor_id='%s' " % icgc_donor_id
			ret2 = search_db(cursor,qry)
			print(ret2)
			# check the full length of the entry
			qry  = "select * from %s " % table
			qry += "where icgc_mutation_id = '%s' " % icgc_mutation_id
			qry += "and icgc_donor_id = '%s' " % icgc_donor_id
			qry += "and icgc_specimen_id = '%s' " % icgc_specimen_id
			qry += "and icgc_sample_id = '%s' " % icgc_sample_id
			ret2 = search_db(cursor,qry)
			[max_depth, max_allele_depth] = [-1,-1]
			[max_id, max_allele_id]  = [-1,-1]
			all_ids = []
			genotype  = []
			path_estimate = []
			for line2 in ret2:
				named_field = dict(list(zip(colnames,line2)))
				all_ids.append( named_field["id"])
				genotype.append(named_field["tumor_genotype"])
				path_estimate.append(named_field["pathogenic_estimate"])
				if named_field["total_read_count"] and max_depth<named_field["total_read_count"]:
					max_depth = named_field["total_read_count"]
					max_id    = named_field["id"]
				# as the second line tiebreaker
				if named_field["mutant_allele_read_count"] and max_allele_depth<named_field["mutant_allele_read_count"]:
					max_allele_depth = named_field["mutant_allele_read_count"]
					max_allele_id    = named_field["id"]
				print(line2)

			if max_id>=0 and max_allele_id>=0:
				other_ids = set(all_ids)
				if max_id>=0: # ids start from 1
					other_ids.remove(max_id)
				elif max_allele_id>=0:
					other_ids.remove(max_allele_id)
				for other_id in other_ids:
					qry = "delete from %s where id = %d" % (table, other_id)
					search_db(cursor,qry)
			elif ct==2 and genotype[0]==genotype[1][::-1]: # hack to reverse a string
				print("tumor genotypes same", genotype)
				qry = "delete from %s where id = %d" % (table, all_ids[1])
				search_db(cursor,qry)
			# here I am losing patience a bit, I guess
			elif set(path_estimate)=={0}:
				print("all estimated irrelevant")
				for id in all_ids:
					qry = "delete from %s where id = %d" % (table, id)
					search_db(cursor,qry)
			else: # I really cannot decide; therefore merge the annotations into the first entry and delete the rest
				qry1 = "update %s set tumor_genotype='%s' where id=%d" % (table, ";".join(genotype), all_ids[0])
				search_db(cursor,qry1)
				for other_id in all_ids[1:]:
					qry2 = "delete from %s where id = %d" % (table, other_id)
					search_db(cursor,qry2)



	cursor.close()
	db.close()


	return



#########################################
if __name__ == '__main__':
	main()
