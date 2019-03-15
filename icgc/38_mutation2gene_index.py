#! /usr/bin/python
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

# mutation2gene is not empty, do this before 37_mutation_to_gene

import time

from icgc_utils.mysql   import  *
from icgc_utils.processes import *
from random import shuffle


#########################################
#########################################
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()
	#########################
	switch_to_db(cursor,"icgc")
	#qry = "create index gene_idx on mutation2gene (gene_symbol)"
	qry = "create index mut_idx on mutation2gene (icgc_mutation_id)"

	search_db(cursor, qry, verbose=True)

	return

#########################################
if __name__ == '__main__':
	main()
