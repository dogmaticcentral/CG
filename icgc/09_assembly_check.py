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
import subprocess
import time

from config import Config
from icgc_utils.common_queries import  *
from icgc_utils.processes import  *
from icgc_utils.annovar import *
from icgc_utils.CrossMap import *


#########################################
#########################################
def main():
	# what is the difference between GRCh37 and hg19 ?
	# https://www.biostars.org/p/123767/
	# The genomic content for the two is identical, except for the mitochondrial contig
	# we do not have  MT, so wither one should be fine
	ref_assemblies = ['hg19', 'GRCh37'] # this is the assembly I would like to see all coords in location tables
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	for table in tables:
		print("======================")
		print(table)
		# get all coords which are not reference
		qry =  "select  distinct assembly from {}". format(table)
		ret = search_db(cursor,qry,verbose=True)
		print(ret)


	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
