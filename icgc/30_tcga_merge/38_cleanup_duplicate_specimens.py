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


from icgc_utils.common_queries   import  *
from config import Config
verbose = True


def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")
	for table in tables:
		# specimens per donor?
		qry  = "select  icgc_donor_id, count(distinct  icgc_specimen_id) ct "
		qry += "from  %s  " % table
		qry += "group by icgc_donor_id having ct>1 order by ct desc"
		ret = search_db(cursor,qry)
		if not ret: continue
		if verbose:
			print("=================================\n",table)
			print("donors_with_multiple_specimens: ", len(ret))
		icgc_donor_ids = [r[0] for r in ret]
		specimen_ids = {}
		for icgc_donor_id in icgc_donor_ids:
			qry  = "select distinct(icgc_specimen_id) "
			qry += "from %s where icgc_donor_id='%s' " % (table,icgc_donor_id)
			ret  = search_db(cursor,qry)
			if ret and len(ret)>1:
				specimen_ids[icgc_donor_id] = [r[0] for r in ret]

		resolve_duplicates(cursor, table, icgc_donor_ids, specimen_ids)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
