#! /usr/bin/python
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
