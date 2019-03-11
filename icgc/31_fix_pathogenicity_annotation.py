#! /usr/bin/python3
# I've managed to screw this one up (missed splice?), so go back and fix

import subprocess
import time, re

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *
from random import shuffle
from config import Config

# this is set literal
mutation_pathogenic = {'missense','frameshift',  'stop_gained', 'inframe',
              'stop_lost', 'inframe_deletion', 'inframe_insertion',
              'start_lost', 'disruptive_inframe_deletion',
               'exon_loss', 'disruptive_inframe_insertion',
              'splice', '5_prime_UTR_premature_start_codon_gain',
              'splice_acceptor', 'splice_region', 'splice_donor'
             }

location_pathogenic = { 'splice', '5_prime_UTR_premature_start_codon_gain',
              'splice_acceptor', 'splice_region', 'splice_donor',
}
#########################################
def fix_pathogenicity(chromosomes, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	for chrom  in chromosomes:
		mutations_table = "mutations_chrom_%s"%chrom
		print()
		print("====================")
		print("processing icgc table ", mutations_table, os.getpid())
		qry = "select icgc_mutation_id, start_position, consequence, pathogenic_estimate from icgc.%s" % mutations_table
		for  icgc_mutation_id, start_position, consequence, p_estimate in search_db(cursor,qry):
			#print icgc_mutation_id, start_position, consequence, p_estimate
			p_estimate_revised = 0
			if consequence:
				for description in mutation_pathogenic:
					if description in consequence:
						p_estimate_revised=1
			if p_estimate_revised == 0: # check if location is splice
				locations_table = "locations_chrom_%s"%chrom
				qry2 = "select transcript_relative from icgc.%s " % locations_table
				qry2 += "where position = %d" % start_position
				# position is the principal key, so there should not be two of those
				ret = search_db(cursor,qry2)
				if not ret:
					search_db(cursor,qry2,verbose=True)
					exit()
				tr_relative = ret[0][0]
				#print "tr_relative", tr_relative
				if tr_relative:
					for description in  location_pathogenic:
						if description in tr_relative:
							p_estimate_revised=1
			#print p_estimate, p_estimate_revised
			if (p_estimate==None) or  (p_estimate !=  p_estimate_revised): # if the revision needed, proceed
				qry3  = "update icgc.%s " %  mutations_table
				qry3 += "set pathogenic_estimate=%d " % p_estimate_revised
				qry3 += "where icgc_mutation_id='%s' " %  icgc_mutation_id
				search_db(cursor,qry3)

	cursor.close()
	db.close()
	return


#########################################
def main():

	chromosomes = [str(i) for i in range(1,23)] + ["X","Y"]
	shuffle(chromosomes)


	number_of_chunks = 8  # myISAM does not deadlock
	parallelize(number_of_chunks, fix_pathogenicity, chromosomes, [])

#########################################
if __name__ == '__main__':
	main()
