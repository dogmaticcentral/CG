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
from icgc_utils.common_queries import *
from config import Config
from random import sample
from math import sqrt

####################################################
def main():


	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	protein_coding, chrom = protein_coding_genes(cursor)

	switch_to_db(cursor, 'icgc')

	print("never mutated")
	count = 0
	for gene in protein_coding:
		qry = "select * from mutation2gene where gene_symbol='%s' limit 1" % gene
		ret = error_intolerant_search(cursor, qry)
		if ret: continue
		print(gene)
		count += 1
	print("count:", count)


	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
