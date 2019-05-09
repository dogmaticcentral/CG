#! /usr/bin/python3
#
# This source code is part of icgc, an ICGC processing pipeline.
# 
# Icgc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Icgc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#
from config import Config
from icgc_utils.common_queries  import  *
from icgc_utils.utils import *
from icgc_utils.processes   import  *

def find_duplicates(cursor, somatic_table, diagnostics=False):

	duplicates = []

	tumor_short = somatic_table.split("_")[0]

	##################
	outf = None
	if diagnostics: outf = open("%s.duplicate_donors.tsv"%tumor_short, "w")
	if outf: outf.write("\t".join(['tumor', 'donor 1', 'donor 2', 'vars 1', 'vars 2', 'common vars', 'pct 1', 'pct2'])+"\n")

	##################
	# we are looking for suspiciously large overlaps in the mutation sets
	# first find the mutation set associated with each mutation id
	qry = "select icgc_donor_id,  group_concat(icgc_mutation_id) from %s " % somatic_table
	qry += "group by icgc_donor_id"
	ret = search_db(cursor, qry, verbose=False)
	if not ret: return None

	mutations = {}
	for icgc_donor_id, mutstr in ret:
		mutations[icgc_donor_id] = set(mutstr.split(","))

	# sanity - the only information that we have available that might eliminate
	# a pair as potential duplicate is the sex of the donor
	sex = {}
	qry = "select icgc_donor_id, donor_sex from %s_donor " % tumor_short
	ret = search_db(cursor, qry, verbose=False)
	if ret: sex = dict(ret)

	# now consider each pair of ids and the overlap between their mutation sets
	donor = list(mutations.keys())
	for i in range(len(mutations)):
		mi =  mutations[donor[i]]
		li = len(mi)
		sex_i = sex[donor[i]] if donor[i] in sex else None
		for j in range(i+1,len(mutations)):
			mj = mutations[donor[j]]
			lj = len(mj)
			sex_j = sex[donor[j]] if donor[j] in sex else None
			if sex_i and sex_j and sex_i.lower()!=sex_j.lower(): continue
			common_vars =  mi.intersection(mj)
			lv = len(common_vars)
			if lv<5: continue
			pcti = float(lv)/li
			pctj = float(lv)/lj
			if pcti<0.4 and pctj<0.4: continue
			#print("\t".join([str(i) for i in [tumor_short, donor[i], donor[j], li, lj, lv]]) )
			if outf:
				outfields = [tumor_short, donor[i], donor[j], li, lj, lv, "%d"%(100*pcti), "%d"%(100*pctj)]
				outf.write("\t".join([str(i) for i in outfields])+"\n" )
			duplicates.append([donor[i], donor[j]])

	if outf: outf.close()

	return duplicates


#########################################
def resolve(cursor, somatic_table, cluster):


	tumor_short = somatic_table.split("_")[0]

	# by this point is should have resolved duplicate specimens, but let's check again
	qry  = "select icgc_donor_id,  group_concat(distinct(icgc_specimen_id)) from %s " % somatic_table
	qry += "where icgc_donor_id in (%s) " % ",".join([quotify(d) for d in cluster])
	qry += "group by icgc_donor_id"
	ret = search_db(cursor, qry, verbose=False)
	specimen_id = dict(ret)
	for icgc_donor_id, spid in specimen_id.items():
		if "," in spid:
			print("{}: multiple specimen ids ({}) for {}".format(somatic_table, spid, icgc_donor_id))
			exit()

	# specimen types?
	qry  = "select icgc_donor_id, specimen_type from %s_specimen " % tumor_short
	qry += "where icgc_specimen_id in (%s)" % ",".join([quotify(d) for d in specimen_id.values()])
	ret  = search_db(cursor, qry)
	specimen_type = dict(ret)

	# total number of variants, depth of sequencing?
	qry  = "select icgc_donor_id, count(*), ifnull(avg(total_read_count),0)  from %s " % somatic_table
	qry += "where icgc_donor_id in (%s) group by icgc_donor_id" % ",".join([quotify(d) for d in cluster])
	ret  = search_db(cursor, qry)

	number_of_variants = dict([r[:2] for r in ret])
	avg_depth = dict([[r[0],r[2]] for r in ret])
	print("\ncluster")
	for donor_id in cluster:
		print(" %10s  %10d  %10.2f  %s" %(donor_id, number_of_variants[donor_id], avg_depth[donor_id], specimen_type[donor_id]))
	print("end cluster")


	primary_tumors = set([donor_id for donor_id in cluster if 'primary' in specimen_type[donor_id].lower()[:len('primary')]])
	metastatic_tumors =  set([donor_id for donor_id in cluster if 'metastatic' in specimen_type[donor_id].lower()[:len('metastatic')]])
	normal_specimens =  set([donor_id for donor_id in cluster if 'normal' in specimen_type[donor_id].lower()[0:len('normal')]])
	other = cluster.difference(primary_tumors).difference(metastatic_tumors).difference(normal_specimens)

	preferred_set = other
	for tumor_type_set in [primary_tumors, metastatic_tumors]:
		if len(tumor_type_set)>0:
			preferred_set = tumor_type_set
			break
	if len(preferred_set)==0: return # not sure how this could happen
	keep_id = None

	if len(preferred_set)==1:
		keep_id = preferred_set.pop()
	else:
		keep_id = sorted(preferred_set,key=lambda x: avg_depth[x], reverse=True)[0]

	cluster.remove(keep_id)
	print("*** keeping", keep_id)
	print("### dropping", cluster)
	remove_string = ",".join([quotify(d) for d in cluster])
	qry = "delete from %s where icgc_donor_id in (%s)" % (somatic_table, remove_string)
	ret  = search_db(cursor, qry)
	if ret:
		search_db(cursor, qry, verbose=True)
		exit()
	print(qry)
	return

#########################################
def cleanup(tables, other_args):
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	# group concat default length is 1024 bytes which is kinda smallish
	qry = "set session group_concat_max_len = 20000000"
	search_db(cursor, qry)

	for somatic_table in tables:
		print("\n==========\n%s"%somatic_table)
		duplicates = find_duplicates(cursor, somatic_table, diagnostics=False)
		if not duplicates or len(duplicates)==0: continue
		clusters_of_duplicates = find_clusters(duplicates)
		for cluster in clusters_of_duplicates:
			resolve(cursor, somatic_table, cluster)
	cursor.close()
	db.close()


#########################################
#########################################
def main():

	#print("disabled ") # comment out to run
	#exit(1)


	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which  somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	table_size = get_table_size(cursor, 'icgc', tables)
	cursor.close()
	db.close()

	tables_sorted = sorted(tables, key=lambda t: table_size[t], reverse=False)
	number_of_chunks = 1
	parallelize(number_of_chunks, cleanup, tables_sorted, [], round_robin=True)


	return



#########################################
if __name__ == '__main__':
	main()
