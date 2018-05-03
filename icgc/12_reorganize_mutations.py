#! /usr/bin/python



# consequence type
# location:
# actual consequence

from icgc_utils.common_queries   import  *

variant_columns = ['icgc_mutation_id', 'chromosome','icgc_donor_id', 'icgc_specimen_id', 'icgc_sample_id',
                   'control_genotype', 'tumor_genotype', 'total_read_count', 'mutant_allele_read_count']

mutation_columns = ['icgc_mutation_id', 'start_position', 'end_position', 'mutation_type',
					'mutated_from_allele', 'mutated_to_allele', 'reference_genome_allele',
					'aa_mutation', 'cds_mutation']

location_columns = ['position', 'gene_relative', 'transcript_relative']

################################################################
# stop_retained: A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
consequence_vocab = ['stop_lost', 'synonymous', 'inframe_deletion', 'inframe_insertion', 'stop_gained',
                     '5_prime_UTR_premature_start_codon_gain',
                     'start_lost', 'frameshift', 'disruptive_inframe_deletion', 'stop_retained',
                     'exon_loss', 'disruptive_inframe_insertion', 'missense']
# this is set literal
pathogenic =  {'stop_lost', 'inframe_deletion', 'stop_gained', '5_prime_UTR_premature_start_codon_gain',
                     'start_lost', 'frameshift', 'disruptive_inframe_deletion',
                     'exon_loss', 'disruptive_inframe_insertion', 'missense'}

# location_vocab[1:4] is gene-relative
# location_vocab[1:4] is transcript-relative
location_vocab = ['intergenic_region', 'intragenic', 'upstream', 'downstream',
                  '5_prime_UTR', 'exon',  'coding_sequence', 'initiator_codon',
                  'splice_acceptor', 'splice_region', 'splice_donor',
                  'intron', '3_prime_UTR', ]

#########################################
def quotify(something):
	if not something:
		return ""
	if type(something)==str:
		return "\'"+something+"\'"
	else:
		return str(something)

#########################################
def reorganize_donor_variants(cursor, table, columns):
	variants_table = table.replace("_temp","")
	donors = get_donors(cursor, table)
	for donor in donors:
		variants  = set([])
		total_entries = 0
		for fields in search_db (cursor, "select * from %s where icgc_donor_id='%s'" % (table, donor)):
			total_entries += 1
			named_field = dict(zip(columns,fields))

			variant_values = []
			for name in variant_columns:
				variant_values.append(quotify(named_field[name]))
			variants.add(",".join(variant_values)) # this is set, getting rid of duplicate info

	for variant in variants:
		qry = "insert into %s (%s) " %(variants_table, ",".join(variant_columns))
		qry += "values (%s) " % variant
		search_db(cursor, qry, verbose=True)

#########################################
def reorganize_mutations(cursor, table, columns):

	mutations = get_mutations(cursor, table)
	for mutation in mutations:

		conseqs   = set([])
		mutation_values = None
		chromosome = None
		for fields in search_db (cursor, "select * from %s where icgc_mutation_id='%s'" % (table, mutation)):

			named_field = dict(zip(columns,fields))

			if not mutation_values: # we need to set this only once
				mutation_values = [quotify(named_field[name]) for name in mutation_columns]
				chromosome = named_field['chromosome']

			# this is not ready to be stored, because we need to work through the consequences
			csq = named_field['consequence_type']
			if csq in consequence_vocab:
				conseqs.add(csq)
			elif csq in location_vocab:
				pass
			elif csq == "":
				pass
			else:
				print "unrecognized consequence field:", csq
				exit()

		if not mutation_values:
			print "mutation values not assigned for %s (!?)" % mutation
			exit()
		if not chromosome:
			print "chromosome not assigned for %s (!?)" % mutation
			exit()

		mutation_values.append(quotify(";".join(list(conseqs))))
		if len(conseqs&pathogenic)>0:
			mutation_values.append("1")
		else:
			mutation_values.append("0")

		# now we are ready to store
		mutation_table = "mutations_chrom_{}".format(chromosome)
		qry = "insert into %s (%s) " %(mutation_table, ",".join(mutation_columns + ['consequence', 'pathogenic_estimate']))
		qry += "values (%s) " % ",".join(mutation_values)
		search_db(cursor, qry, verbose=True)
		exit()

#########################################
def reorganize_locations(cursor, table, columns):

	chromosomes = [str(i) for i in range(1,23)] + ["X", "Y", "MT"]

	for chromosome in chromosomes:
		qry =  "select distinct start_position from %s  where chromosome='%s'" % (table, chromosome)
		positions = [ret[0] for ret in  search_db (cursor, qry)]
		for position in positions:
			gene_relative       = set([])
			transcript_relative = set([])
			qry =  "select * from %s where chromosome='%s' and start_position=%d" % (table, chromosome, position)
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

			# now we are ready to store
			location_values = [str(position), quotify(";".join(gene_relative)), quotify(";".join(transcript_relative))]
			location_table = "locations_chrom_{}".format(chromosome)
			qry = "insert into %s (%s) " %(location_table, ",".join(location_columns))
			qry += "values (%s) " % ",".join(location_values)
			print qry

#########################################
def reorganize(cursor, table, columns):
	# line by line: move id info into new table
	# for mutation and location check if the info exists; if not make new entry
	print "===================="
	print "reorganizing ", table
	###################################
	#reorganize_donor_variants(cursor, table, columns)
	reorganize_mutations(cursor, table, columns)
	#reorganize_locations(cursor, table, columns)


	exit()



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

	# the tables should all have the same columns
	qry = "select column_name from information_schema.columns where table_name='%s'"%tables[0]
	columns = [field[0] for field in  search_db(cursor,qry)]

	switch_to_db(cursor,"icgc")
	# enable if run for the first time
	#sanity_checks(cursor, tables)

	tables = ["EOPC_simple_somatic_temp"]

	for table in tables:
		reorganize (cursor, table, columns)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
