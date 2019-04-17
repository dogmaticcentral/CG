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

from icgc_utils.mysql import *
from config import Config

#########################################
#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_donor'"

	tables = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor,'icgc')
	for table in tables:
		print("======================")
		print(table)
		tumor_short = table.split("_")[0]
		variants_table = table.replace("donor","simple_somatic_temp")
		# get all coords which are not reference
		qry =  "select  icgc_donor_id, submitted_donor_id from {} ". format(table)
		qry += "where icgc_donor_id not like 'DOT%'"
		ret = search_db(cursor,qry,verbose=True)
		number_of_ids = len(ret) if ret else 0
		print("number of donor ids:", number_of_ids)
		if number_of_ids==0: continue
		donors_with_no_variants = 0
		for icgc_donor_id, submitted_donor_id in ret:
			qry = "select count(1) from %s where icgc_donor_id='%s'" % (variants_table, icgc_donor_id)
			ret2 = search_db(cursor,qry,verbose=True)
			number_of_variants = ret2[0][0]
			#print(icgc_donor_id, submitted_donor_id, "number_of_variants: ", number_of_variants)
			if number_of_variants==0: donors_with_no_variants +=1
		print("icgc donor ids with no variants:",donors_with_no_variants )
	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
