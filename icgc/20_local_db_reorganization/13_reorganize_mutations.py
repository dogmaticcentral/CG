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

from config import Config
from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *


#  'aa_mutation',  'consequence_type', and 'pathogenic_estimate'  will be filled separately
mutation_columns = ['icgc_mutation_id', 'start_position', 'end_position', 'assembly',
					'mutation_type', 'mutated_from_allele', 'mutated_to_allele', 'reference_genome_allele']


################################################################
# stop_retained: A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
consequence_vocab = ['stop_lost', 'synonymous', 'inframe_deletion', 'inframe_insertion', 'stop_gained',
                     '5_prime_UTR_premature_start_codon_gain',
                     'start_lost', 'start_gained', 'frameshift', 'disruptive_inframe_deletion', 'stop_retained',
                     'exon_loss', 'disruptive_inframe_insertion', 'missense']

# location_vocab[1:4] is gene-relative
# location_vocab[1:4] is transcript-relative
location_vocab = ['intergenic_region', 'intragenic', 'upstream', 'downstream',
                  '5_prime_UTR', 'exon',  'coding_sequence', 'initiator_codon',
                  'splice_acceptor', 'splice_region', 'splice_donor',
                  'intron', '3_prime_UTR', ]

# this is set literal
pathogenic = {'stop_lost', 'inframe_deletion', 'inframe_insertion', 'stop_gained',
                    '5_prime_UTR_premature_start_codon_gain',
                     'start_lost', 'start_gained', 'frameshift', 'disruptive_inframe_deletion',
                     'exon_loss', 'disruptive_inframe_insertion', 'missense',
                     'splice_acceptor', 'splice_region', 'splice_donor',
                      'inframe'   # there is no way we can know at this level whether an inframe change is nondsisruptive
									# more likely it is than not
             }


#########################################
def insert (cursor, table, columns, values):

	nonempty_values = []
	corresponding_columns = []
	for i in range(len(values)):
		if not values[i] or  values[i] == "": continue
		nonempty_values.append(values[i])
		corresponding_columns.append(columns[i])
	qry = "insert into %s (%s) " %(table, ",".join(corresponding_columns))
	qry += "values (%s) " % ",".join(nonempty_values)
	search_db(cursor, qry)


#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./icgc_utils/kernprof.py -l 12_reorganize_variants.py
# followed by
# python -m line_profiler 12_reorganize_variants.py.lprof
# see here https://github.com/rkern/line_profiler#line-profiler
# the reason I am using local kernprof.py is that I don't know where pip
# installed its version (if anywhere)
# @profile
def reorganize_mutations(cursor, table, columns):
	# reorganize = divide into three tables: variants(per user), mutations, and locations
	mutations = get_mutations(cursor, table)
	totmut = len(mutations)
	print("\t\t\t total mutations:", totmut)
	ct = 0
	time0 = time.time()
	for mutation in mutations:
		ct += 1
		if ct%10000 == 0:
			print("\t\t\t %10s  %6d  %d%%  %ds" % (table, ct, float(ct)/totmut*100, time.time()-time0))
			time0 = time.time()
		mutation_already_seen = False
		conseqs   = set([])
		aa_mutations = set([])
		mutation_values = None
		chromosome = None
		mutation_table = None

		# this hinges on index on the *simple_somatic_temp
		# qry  = "create index mut_gene_idx on %s (icgc_mutation_id, gene_affected)" % simple_somatic_temp_table
		qry  = "select * from %s where icgc_mutation_id='%s' " % (table, mutation)
		qry += "and gene_affected is not null and gene_affected !='' "
		ret  = search_db (cursor, qry)

		if not ret: continue
		for fields in ret:

			named_field = dict(list(zip(columns,fields)))

			if not mutation_values: # we need to set this only once
				mutation_values = [quotify(named_field[name]) for name in mutation_columns]
				chromosome = named_field['chromosome']
				mutation_table = "mutations_chrom_{}".format(chromosome)
				if entry_exists(cursor, "icgc", mutation_table, "icgc_mutation_id", quotify(mutation)):
					mutation_already_seen = True
					continue

			# aa_mutation
			aa = named_field['aa_mutation']
			if aa and  aa!="":
				transcript =  named_field['transcript_affected']
				if not transcript: transcript="unk"
				aa_mutations.add("{}:{}".format(transcript,aa))

			# consequences
			csq = named_field['consequence_type']
			if csq in consequence_vocab:
				conseqs.add(csq)
			elif csq in location_vocab:
				if "splice" in csq.lower():
					conseqs.add(csq)
			elif csq == "":
				pass
			else:
				print("unrecognized consequence field:", csq)
				exit()

		if mutation_already_seen: continue

		if not mutation_values:
			print("mutation values not assigned for %s (!?)" % mutation)
			exit()
		if not chromosome:
			print("chromosome not assigned for %s (!?)" % mutation)
			exit()

		mutation_values.append(quotify(";".join(list(aa_mutations))))
		mutation_values.append(quotify(";".join(list(conseqs))))
		if len(conseqs&pathogenic)>0:
			mutation_values.append("1")
		else:
			mutation_values.append("0")

		# now we are ready to store
		insert(cursor, mutation_table, mutation_columns + ['aa_mutation','consequence', 'pathogenic_estimate'], mutation_values)


#########################################
def reorganize(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:

		# the tables should all have the same columns
		qry = "select column_name from information_schema.columns where table_name='%s'"%table
		columns = [field[0] for field in  search_db(cursor,qry)]
		# line by line: move id info into new table
		# for mutation and location check if the info exists; if not make new entry
		time0 = time.time()
		print("====================")
		print("reorganizing mutations from", table, os.getpid())
		reorganize_mutations(cursor, table, columns)
		time1 = time.time()
		print(("\t\t done in %.3f mins" % (float(time1-time0)/60)), os.getpid())

	cursor.close()
	db.close()

	return


#########################################
#########################################
def main():

	print("disabled")
	exit()

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	number_of_chunks = 8  # myISAM does not deadlock
	parallelize(number_of_chunks, reorganize, tables, [])



#########################################
if __name__ == '__main__':
	main()
