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
from icgc_utils.mysql import *
from config import Config

####################################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, 'icgc')
	for pthwy_id, name in hard_landing_search(cursor, "select * from reactome_pathways"):
		if "'" in name:
			print("*{}*".format(name))
			continue
		name = name.strip()
		qry = "update reactome_pathways set name='%s' where reactome_pathway_id='%s'" % (name, pthwy_id)
		error_intolerant_search(cursor, qry)

	cursor.close()
	db.close()

	return


#########################################
if __name__ == '__main__':
	main()
