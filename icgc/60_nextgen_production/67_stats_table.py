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
import time

from icgc_utils.common_queries import *
from icgc_utils.processes import *
from config import Config
from random import sample
from math import sqrt

def make_tables(cursor, db_name):
	switch_to_db(cursor, db_name)

	table = 'stats_description'
	if not check_table_exists(cursor, db_name, table):
		qry = ""
		qry += "  CREATE TABLE  %s (" % table
		qry += "  	 stats_id VARCHAR (20) NOT NULL, "
		qry += "     description TEXT NOT NULL, "
		qry += "     parameters TEXT NOT NULL, "
		qry += "     stats TEXT NOT NULL, "
		qry += "	 PRIMARY KEY (stats_id) "
		qry += ") ENGINE=MyISAM"
		search_db(cursor, qry, verbose=True)

	table = 'stats'
	if not check_table_exists(cursor, db_name, table):
		qry = ""
		qry += "  CREATE TABLE  %s (" % table
		qry += "     id INT NOT NULL AUTO_INCREMENT, "
		qry += "  	 stats_id VARCHAR (20) NOT NULL, "
		qry += "  	 parameters TEXT NOT NULL, "
		qry += "  	 stats TEXT  NOT NULL, "
		qry += "	 PRIMARY KEY (id) "
		qry += ") ENGINE=MyISAM"
		search_db(cursor, qry, verbose=True)

#########################################
#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	db_name =  "icgc"
	make_tables(cursor, db_name)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
