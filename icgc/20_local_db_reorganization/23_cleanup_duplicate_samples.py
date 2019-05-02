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

def count_entries(cursor, somatic_table, icgc_sample_id):
	qry = "select count(*) from {} where icgc_sample_id='{}' ".format(somatic_table, icgc_sample_id)
	ret = search_db(cursor,qry)
	if not ret or ret[0][0]==0: return 0
	return ret[0][0]

def find_sample_id_with_max_entries(sample_ids, entries_per_sample):
	max_count = 0
	max_sample_id = None
	for sample_id in sample_ids:
		if not sample_id in entries_per_sample: continue
		if max_count<entries_per_sample[sample_id]:
			max_count=entries_per_sample[sample_id]
			max_sample_id = sample_id
	return max_sample_id

#########################################
#########################################
def main():

	print("disabled - change to merge rather than delete - the greatest depth or some such") # comment out to run
	exit(1)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which  somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,"icgc")
	for somatic_table in tables:
		icgc_donor_ids = [r[0] for r in search_db(cursor, "select distinct icgc_donor_id from %s"%somatic_table)]
		problematic = []
		sample_ids = {}
		for icgc_donor_id in icgc_donor_ids:
			qry  = "select distinct(icgc_sample_id) "
			qry += "from %s where icgc_donor_id='%s' " % (somatic_table,icgc_donor_id)
			ret  = search_db(cursor,qry)
			if ret and len(ret)>1:
				problematic.append(icgc_donor_id)
				sample_ids[icgc_donor_id] = [r[0] for r in ret]
		if len(problematic)==0:
			print("%s has no duplicate samples in the variant table"% somatic_table)
			continue

		print("%s has %d donor ids with duplicate samples ids in the variant table" % (somatic_table,len(problematic)))

		# most of these are innocuous, with normal sample not appearing in the variants table
		# this, however is not always the case
		for icgc_donor_id in problematic:
			qry  = "select icgc_sample_id, count(*) as c  from %s " % somatic_table
			qry += "where icgc_donor_id = '%s' and reliability_estimate=1 " % icgc_donor_id
			qry += "group by  icgc_sample_id"
			ret2 = search_db(cursor,qry)
			entries_per_sample = dict(ret2)
			removable_ids = sample_ids[icgc_donor_id]
			keep_id = find_sample_id_with_max_entries(removable_ids, entries_per_sample)
			if keep_id:
				removable_ids.remove(keep_id)
				removable_ids_string = ",".join(["'%s'"%rid for rid in removable_ids])
				qry = "delete from {} where icgc_sample_id in ({})".format(somatic_table, removable_ids_string)
				search_db(cursor,qry)


	cursor.close()
	db.close()


	return



#########################################
if __name__ == '__main__':
	main()
