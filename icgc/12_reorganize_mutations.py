#! /usr/bin/python


import time

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *

variant_columns = ['icgc_mutation_id', 'chromosome','icgc_donor_id', 'icgc_specimen_id', 'icgc_sample_id',
                   'submitted_sample_id','control_genotype', 'tumor_genotype', 'total_read_count', 'mutant_allele_read_count']

# we'll take care of 'aa_mutation' and 'consequence_type will be handled separately
mutation_columns = ['icgc_mutation_id', 'start_position', 'end_position', 'mutation_type',
					'mutated_from_allele', 'mutated_to_allele', 'reference_genome_allele']


location_columns = ['position', 'gene_relative', 'transcript_relative']

################################################################
# stop_retained: A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
consequence_vocab = ['stop_lost', 'synonymous', 'inframe_deletion', 'inframe_insertion', 'stop_gained',
                     '5_prime_UTR_premature_start_codon_gain',
                     'start_lost', 'frameshift', 'disruptive_inframe_deletion', 'stop_retained',
                     'exon_loss', 'disruptive_inframe_insertion', 'missense']

# location_vocab[1:4] is gene-relative
# location_vocab[1:4] is transcript-relative
location_vocab = ['intergenic_region', 'intragenic', 'upstream', 'downstream',
                  '5_prime_UTR', 'exon',  'coding_sequence', 'initiator_codon',
                  'splice_acceptor', 'splice_region', 'splice_donor',
                  'intron', '3_prime_UTR', ]

# this is set literal
pathogenic = {'stop_lost', 'inframe_deletion', 'stop_gained', '5_prime_UTR_premature_start_codon_gain',
                     'start_lost', 'frameshift', 'disruptive_inframe_deletion',
                     'exon_loss', 'disruptive_inframe_insertion', 'missense',
                     'splice_acceptor', 'splice_region', 'splice_donor'
             }

#########################################
def quotify(something):
	if not something:
		return ""
	if type(something)==str:
		return "\'"+something+"\'"
	else:
		return str(something)

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
def reorganize_donor_variants(cursor, table, columns):

	variants_table = table.replace("_temp","")
	donors = get_donors(cursor, table)
	for donor in donors:
		variants  = set([])
		total_entries = 0
		qry  = "select * from %s where icgc_donor_id='%s' " % (table, donor)
		qry += "and gene_affected is not null and gene_affected !='' "
		ret  = search_db (cursor, qry)
		if not ret: continue # it happens, check "select * from ALL_simple_somatic_temp where icgc_donor_id='DO282'"
		#	print search_db (cursor, qry, verbose=True)
		#	exit()
		for fields in ret:
			total_entries += 1
			named_field = dict(zip(columns,fields))
			variant_values = []
			for name in variant_columns:
				variant_values.append(quotify(named_field[name]))
			variants.add(",".join(variant_values)) # set => getting rid of duplicates

		for variant in variants:
			insert(cursor, variants_table, variant_columns, variant.split(","))

#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./icgc_utils/kernprof.py -l 12_reorganize_mutations.py
# followed by
# python -m line_profiler 12_reorganize_mutations.py.lprof
# see here https://github.com/rkern/line_profiler#line-profiler
# the reason I am using local kernprof.py is that I don't know where pip
# installed its version (if anywhere)
# @profile
def reorganize_mutations(cursor, table, columns):
	# reorganize = divide into three tables: variants(per user), mutations, and locations
	mutations = get_mutations(cursor, table)
	totmut = len(mutations)
	print "\t\t\t total mutations:", totmut
	ct = 0
	time0 = time.time()
	for mutation in mutations:
		ct += 1
		if ct%10000 == 0:
			print "\t\t\t %10s  %6d  %d%%  %ds" % (table, ct, float(ct)/totmut*100, time.time()-time0)
			time0 = time.time()
		mutation_already_seen = False
		conseqs   = set([])
		aa_mutations = set([])
		mutation_values = None
		chromosome = None
		mutation_table = None

		# this hinges on index on the *simple_somatic_temp
		# qry  = "create index mut_gene_idx on %s (icgc_mutation_id, gene_affected)" % mutations_table
		qry  = "select * from %s where icgc_mutation_id='%s' " % (table, mutation)
		qry += "and gene_affected is not null and gene_affected !='' "
		ret  = search_db (cursor, qry)

		if not ret: continue
		for fields in ret:

			named_field = dict(zip(columns,fields))

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
				print "unrecognized consequence field:", csq
				exit()

		if mutation_already_seen: continue

		if not mutation_values:
			print "mutation values not assigned for %s (!?)" % mutation
			exit()
		if not chromosome:
			print "chromosome not assigned for %s (!?)" % mutation
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
def reorganize_locations(cursor, table, columns):

	chromosomes = [str(i) for i in range(1,23)] + ["X", "Y", "MT"]

	for chromosome in chromosomes:
		location_table = "locations_chrom_{}".format(chromosome)
		qry =  "select distinct start_position from %s  where chromosome='%s' " % (table, chromosome)
		qry += "and gene_affected is not null and gene_affected !='' "
		ret  = search_db (cursor, qry)
		if not ret: continue
		positions = [r[0] for r in ret]
		for position in positions:

			if entry_exists(cursor, "icgc", location_table, "position", position): continue

			gene_relative       = set([])
			transcript_relative = set([])
			# this hinges on
			# qry  = "create index chrom_pos_idx on %s (chromosome, start_position)" % mutations_table
			qry =  "select * from %s where chromosome='%s' and start_position=%d " % (table, chromosome, position)
			qry += "and gene_affected is not null and gene_affected !='' "
			for fields in search_db (cursor,qry):
				named_field = dict(zip(columns,fields))
				# this is not ready to be stored, because we need to work through the consequences
				gene   = named_field['gene_affected']
				tscrpt = named_field['transcript_affected']
				csq = named_field['consequence_type']
				if csq in consequence_vocab:
					pass
				elif csq == location_vocab[0]: # intergenic
					pass
				elif csq in location_vocab[1:4]: # gene-relative
					gene_relative.add("{}:{}".format(gene,csq))
				elif csq in location_vocab[4:]: # transcript-relative
					gene_relative.add("{}:{}".format(gene,"intragenic"))
					transcript_relative.add("{}:{}".format(tscrpt,csq))
				elif csq == "":
					pass
				else:
					print "unrecognized consequence field:", csq
					exit()
				# gene and transcript are listed, but no relative position is specified:
				if len(gene_relative)==0 and gene and 'ENSG' in gene:
					gene_relative.add("{}:{}".format(gene,"unk"))
				if len(transcript_relative)==0 and tscrpt and 'ENST' in tscrpt:
					transcript_relative.add("{}:{}".format(tscrpt,"unk"))

			# now we are ready to store
			location_values = [str(position), quotify(";".join(gene_relative)), quotify(";".join(transcript_relative))]
			insert (cursor, location_table, location_columns, location_values)

#########################################
def reorganize(tables, other_args):

	print "disabled"
	exit()


	db     = connect_to_mysql()
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:

		# the tables should all have the same columns
		qry = "select column_name from information_schema.columns where table_name='%s'"%table
		columns = [field[0] for field in  search_db(cursor,qry)]
		# line by line: move id info into new table
		# for mutation and location check if the info exists; if not make new entry
		time0 = time.time()
		print "===================="
		print "reorganizing ", table, os.getpid()
		###############
		print "\t variants", os.getpid()
		#reorganize_donor_variants(cursor, table, columns)
		time1 = time.time()
		print ("\t\t done in %.3f mins" % (float(time1-time0)/60)), os.getpid()
		###############
		print "\t mutations", os.getpid()
		reorganize_mutations(cursor, table, columns)
		time2 = time.time()
		print ("\t\t done in %.3f mins" % (float(time2-time1)/60)), os.getpid()
		###############
		print "\t locations", os.getpid()
		#reorganize_locations(cursor, table, columns)
		time3 = time.time()
		print ("\t\t done in %.3f mins" % (float(time3-time2)/60)), os.getpid()

		print ("\t overall time for %s: %.3f mins" % (table, float(time3-time0)/60)), os.getpid()


	cursor.close()
	db.close()

	return


#########################################
#########################################
def main():

	db     = connect_to_mysql()
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	#tables  = ["BRCA_simple_somatic_temp"]
	#number_of_chunks = 1
	number_of_chunks = 8  # myISAM does not deadlock
	parallelize(number_of_chunks, reorganize, tables, [])





#########################################
if __name__ == '__main__':
	main()
