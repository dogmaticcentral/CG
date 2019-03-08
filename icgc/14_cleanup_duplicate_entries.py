#! /usr/bin/python3
from config import Config
from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *


#########################################
#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./icgc_utils/kernprof.py -l 14_cleanup_duplicate_entries.py
# followed by
# python3 -m line_profiler 14_cleanup_duplicate_entries.py.lprof
# see here https://github.com/rkern/line_profiler#line-profiler
# the reason I am using local kernprof.py is that I don't know where pip
# installed its version (if anywhere)
#@profile
def remove_duplicates(table_rows, other_args):

	table  = other_args[0]
	colnames = other_args[1]
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")

	# loop over all duplicate entries
	for line in table_rows:
		[mega_id, ct] = line
		[icgc_mutation_id, icgc_donor_id, icgc_specimen_id, icgc_sample_id] = mega_id.split("_")

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
			# first see what is the best total read count that we have
			if named_field["total_read_count"] and max_depth<named_field["total_read_count"]:
				max_depth = named_field["total_read_count"]
				max_id    = named_field["id"]
			# as the second line tiebreaker -- this could only happen if the total read count is null
			if named_field["mutant_allele_read_count"] and max_allele_depth<named_field["mutant_allele_read_count"]:
				max_allele_depth = named_field["mutant_allele_read_count"]
				max_allele_id    = named_field["id"]

		if max_id>=0 or max_allele_id>=0:
			other_ids = set(all_ids)
			if max_id>=0: # ids start from 1
				other_ids.remove(max_id)
			#print("max depth %d found at %d, removing"%( max_depth, max_id),other_ids )
			elif max_allele_id>=0:
				other_ids.remove(max_allele_id)
			#print("max allele depth %d found at %d, removing"%(max_allele_depth, max_allele_id),other_ids )
			qry = "delete from %s where id in (%s)" % (table, ",".join([str(other_id) for other_id in other_ids]))
			search_db(cursor,qry)

		# we do not have the info about the depth of the sequencing
		elif ct==2 and genotype[0]==genotype[1][::-1]: # hack to reverse a string
			print("tumor genotypes same", genotype)
			qry = "delete from %s where id = %d" % (table, all_ids[1])
			search_db(cursor,qry)

		# here I am losing patience a bit, I guess
		elif set(path_estimate)=={0}:
			print("all entries estimated irrelevant (in terms of pathogenicity)")
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

	print ("\tprocess {} exiting; worked on table {}".format(get_process_id(), table))

	return

#########################################
#########################################
def main():

	#print("disabled")
	#exit()

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,"icgc")
	for table in tables:
		print("\n====================")
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
			print("\tno duplicates found in", table)
			continue
		print("\t%s has %d duplicates" % (table, len(ret)))
		number_of_chunks = 20  # myISAM does not deadlock
		processes = parallelize(number_of_chunks, remove_duplicates, ret, [table, colnames])
		wait_join(processes)

 
	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
