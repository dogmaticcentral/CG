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


#########################################
#########################################
# produce table of the format
# tumor short | tumor long | number of patients | avg number of mutations per patient |
#  number of patients with mutated rpl5 (%of patients; number of genes which are seen mutated in the same or bigger number of patients)
#  | ditto for rp111

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

	total_donors = 0
	donors_with_multiple_specimens = 0
	donors_with_multiple_samples = 0
	for table in tables:
		tumor_short = table.split("_")[0]
		fields = [tumor_short]
		if verbose: print("\n=================================")
		if verbose: print(table)

		# total number of donors?
		qry    = "select count(distinct icgc_donor_id) from  %s " % table
		donors = search_db(cursor,qry)[0][0]
		total_donors += donors
		if verbose: print("total donors: ", donors)

		# specimens per donor?
		qry  = "select  icgc_donor_id, count(distinct  icgc_specimen_id) ct "
		qry += "from  %s  " % table
		qry += "group by icgc_donor_id having ct>1 order by ct desc"
		ret = search_db(cursor,qry)
		if ret:
			donors_with_multiple_specimens += len(ret)
			if verbose: print("donors_with_multiple_specimens: ", len(ret))
		else:
			if verbose: print("donors_with_multiple_specimens: 0")

		# samples per donor?
		qry  = "select  icgc_donor_id, count(distinct  icgc_sample_id) ct "
		qry += "from  %s  " % table
		qry += "group by icgc_donor_id having ct>1 order by ct desc"
		ret = search_db(cursor,qry)
		if ret:
			donors_with_multiple_samples += len(ret)
			if verbose: print("donors_with_multiple_samples: ", len(ret))
		else:
			if verbose: print("donors_with_multiple_samples: 0")

		#total varinats?
		qry  = "select count(*) from %s" % table
		variants = search_db(cursor,qry)[0][0]
		if verbose: print("variants: ", variants)


	print("total_donors:", total_donors)
	print("donors_with_multiple_specimen labels:", donors_with_multiple_specimens)
	print("donors_with_multiple_sample labels:", donors_with_multiple_samples)
	print()


	for chrom in [str(i) for i in range(1,23)] + ['X','Y']:
		table = "mutations_chrom_%s" % chrom
		qry = "select count(*) from %s where consequence like '%%missense%%'" % table
		miss = search_db(cursor,qry)[0][0]
		qry = "select count(*) from %s where consequence like '%%frameshift%%'" % table
		frm  = search_db(cursor,qry)[0][0]

		print("chromosome %2s   missense: %6d  frameshift: %6d   ratio: %.2f " % (chrom, miss, frm, float(frm)/miss))



	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
