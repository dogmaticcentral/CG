#!/usr/bin/python
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
from icgc_utils.mysql   import  *
from icgc_utils.icgc   import  *


#########################################
#########################################
def main():

	#print("disabled")
	#exit()

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")

	# indices on simple somatic temp
	qry  = "show tables like 'view%'"
	views = [field[0] for field in  search_db(cursor,qry)]
	for view in views:
		print(view)
		qry = "drop view " + view
		search_db(cursor, qry, verbose=True)

	qry  = "show tables like 'temp%'"
	temps = [field[0] for field in  search_db(cursor,qry)]
	for temp in temps:
		print(temp)
		qry = "drop table " + temp
		search_db(cursor, qry, verbose=True)
		
	

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()

