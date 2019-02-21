
from icgc_utils.mysql   import  *

#########################################
def gnomad_mutations (cursor, gene_symbol):

	mutations = []

	chromosome = find_chromosome(cursor, gene_symbol)
	#column_names
	colnames = get_column_names(cursor,"gnomad","gnomad_freqs_chr_1")

	# brute force approach seems to be fast enough for a single gene
	qry = "select * from gnomad.gnomad_freqs_chr_{} where consequences like '%|{}|%' ".format(chromosome, gene_symbol)
	ret = search_db(cursor,qry)
	if not ret:
		print("nothing found for {}, chromosome {}".format(gene_symbol, chromosome))
		exit()
	for line in ret:
		named_fields = dict(list(zip(colnames,line)))
		relevant_variants = []
		for description in named_fields['consequences'].split(","):
			if not 'RPL5' in description: continue
			if not 'missense' in description: continue
			description_field = description.split("|")
			# I don't have ensembl info here - in a more through implementation one should
			# at least go for annotator here
			# for now just hope that the uniprot is canonical
			if len(description_field[2])==0: continue
			relevant_variants.append(description)
		if len(relevant_variants)==0: continue
		if float(named_fields['variant_count'])<2: continue
		freqency = float(named_fields['variant_count'])/named_fields['total_count']
		#print "%.1e" % freqency,
		for description in relevant_variants:
			description_field = description.split("|")
			#print "  ", description_field[7], description_field[8], # example:  280 V/A
			mutations.append(description_field[8].split("/")[0] + description_field[7] + description_field[8].split("/")[1])
		#print

	return list(set(mutations))


#########################################
def transcript_location_cleanup(cursor, loc, gene_stable_id):
	if not loc: return ""
	if loc== "": return loc
	location = {};
	for enst_loc in loc.split(";"):
		[e, c] = enst_loc.split(":")
		location[e] = c
	enst_canonical = list_of_transcript_ids_2_canonical_transcript_id(cursor, list(location.keys()), gene_stable_id)
	if not enst_canonical: return loc
	return location[enst_canonical]

#########################################
def aa_change_cleanup(cursor, aa_change):
	if not aa_change: return ""
	if aa_change=="": return aa_change
	change = {};
	for enst_change in aa_change.split(";"):
		[e, c] = enst_change.split(":")
		change[e] = c
	enst_canonical = list_of_transcript_ids_2_canonical_transcript_id(cursor, list(change.keys()))
	if not enst_canonical or not enst_canonical in change: return aa_change
	return change[enst_canonical]

#########################################
def find_53_status(cursor, tumor_short, specimen):
	# g = gene
	# m = mutation
	# v = variant
	qry  = "select g.icgc_mutation_id, v.pathogenic_estimate, m.consequence, m.aa_mutation, l.transcript_relative "
	qry += "from mutation2gene g, %s_simple_somatic v, mutations_chrom_17 m , locations_chrom_17 l " % (tumor_short)
	qry += "where g.gene_symbol='TP53' "
	qry += "and v.icgc_specimen_id = '%s' "  % specimen
	qry += "and g.icgc_mutation_id = v.icgc_mutation_id "
	qry += "and g.icgc_mutation_id = m.icgc_mutation_id "
	qry += "and m.start_position = l.position "
	ret = search_db(cursor,qry)
	if not ret: return ["wt",""]

	impact_estimate =  "benign"
	cons = []
	for line in ret:
		if line[1]==1:  impact_estimate = "pathogenic"
		if line[2]==None:
			# canonical transcript is ENST00000269305
			if line[-1]!=None: cons.append(transcript_location_cleanup(cursor,line[-1],'ENSG00000141510'))
			continue
		aa_change = aa_change_cleanup(cursor, line[3])
		if aa_change and aa_change != "":
			cons.append("%s:%s"%(line[2], aa_change))
		else:
			cons.append(line[2])
	return [impact_estimate, ";".join(cons)]


#########################################
def protein_coding_genes(cursor):
	standard_chromosomes = [str(i) for i in range(23)] + ['X','Y']
	genes = []
	chrom = {}
	qry  = "select approved_symbol, chromosome from icgc.hgnc "
	qry += "where locus_group='protein-coding gene'"
	for gene,chr in search_db(cursor,qry):
		if not chr in standard_chromosomes: continue
		genes.append(gene)
		chrom[gene] = chr

	return genes, chrom

#########################################
def co_ocurrence_raw(cursor, somatic_table, gene1, gene2):

	qry =  "select g1.gene_symbol,  g2.gene_symbol, s1.icgc_donor_id, s1.submitted_sample_id  "
	qry += "from mutation2gene g1, mutation2gene g2,  %s s1,  %s s2  " % (somatic_table, somatic_table)
	qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
	qry += "and s1.icgc_mutation_id=g1.icgc_mutation_id and g1.gene_symbol='%s' " % gene1
	qry += "and s2.icgc_mutation_id=g2.icgc_mutation_id and g2.gene_symbol='%s' " % gene2
	qry += "and s1.pathogenic_estimate=1 and s1.reliability_estimate=1  "
	qry += "and s2.pathogenic_estimate=1 and s2.reliability_estimate=1 "
	return search_db(cursor,qry)

#########################################
def quotify(something):
	if not something:
		return ""
	if type(something)==str:
		return "\'"+something+"\'"
	else:
		return str(something)


#########################################
def co_ocurrence_w_group_count(cursor, somatic_table, gene1, other_genes):
	qry =  "select count(distinct s1.icgc_donor_id) ct "
	qry += "from mutation2gene g1, mutation2gene g2,  %s s1,  %s s2  " % (somatic_table, somatic_table)
	qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
	qry += "and s1.icgc_mutation_id=g1.icgc_mutation_id and g1.gene_symbol='%s' " % gene1
	group_string = (",".join([quotify(gene2) for gene2 in other_genes]))
	qry += "and s2.icgc_mutation_id=g2.icgc_mutation_id and g2.gene_symbol in (%s) " % group_string
	qry += "and s1.pathogenic_estimate=1 and s1.reliability_estimate=1 "
	qry += "and s2.pathogenic_estimate=1 and s2.reliability_estimate=1 "

	ret = search_db(cursor,qry)

	if not ret:
		search_db(cursor,qry,verbose=True)
		exit()
	return ret[0][0]


#########################################
def co_ocurrence_count(cursor, somatic_table, gene1, gene2):

	if True: # this is still  twice as fast as the search below
		# are we running thruoght the same row twice bcs of s1 <-> s2?
		# still, distinct should get rid of double counting
		qry =  "select count(distinct s1.icgc_donor_id) ct "
		qry += "from mutation2gene g1, mutation2gene g2,  %s s1,  %s s2  " % (somatic_table, somatic_table)
		qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
		qry += "and s1.icgc_mutation_id=g1.icgc_mutation_id and g1.gene_symbol='%s' " % gene1
		qry += "and s2.icgc_mutation_id=g2.icgc_mutation_id and g2.gene_symbol='%s' " % gene2
		qry += "and s1.pathogenic_estimate=1 and s1.reliability_estimate=1  "
		qry += "and s2.pathogenic_estimate=1 and s2.reliability_estimate=1 "
		ret = search_db(cursor,qry)
		if not ret:
			search_db(cursor,qry,verbose=True)
			exit()
		return ret[0][0]
	else:
		donors = {}
		for gene in [gene1, gene2]:
			qry  = "select distinct s.icgc_donor_id "
			qry += "from mutation2gene g,  %s s  " % (somatic_table)
			qry += "where s.icgc_mutation_id=g.icgc_mutation_id and g.gene_symbol='%s' " % gene
			qry += "and s.pathogenic_estimate=1 and s.reliability_estimate=1  "
			ret  = search_db(cursor,qry)
			if not ret:
				search_db(cursor,qry,verbose=True)
				exit()
			donors[gene] = [r[0] for r in ret]

		return len(set(donors[gene1])&set(donors[gene2]))

#########################################
def patients_per_gene_breakdown(cursor, table):

	# this hinges on s.icgc_mutation_id=g.icgc_mutation_id
	# having icgc_mutation_id indexed both on s and g:
	qry  = "select g.gene_symbol symbol, count(distinct  s.icgc_donor_id) ct "
	qry += "from mutation2gene g, %s s  " % table
	qry += "where s.icgc_mutation_id=g.icgc_mutation_id and s.pathogenic_estimate=1  "
	qry += "and s.reliability_estimate=1 "
	qry += "group by symbol"
	ret = search_db(cursor,qry)
	if not ret:
		search_db(cursor,qry, verbose=True)
		exit()
	return dict(ret)

#########################################
def patients_with_muts_in_gene_group(cursor, table, gene_list):

	# this hinges on s.icgc_mutation_id=g.icgc_mutation_id
	# having icgc_mutation_id indexed both on s and g:
	qry  = "select count(distinct  s.icgc_donor_id) ct "
	qry += "from mutation2gene g, %s s  " % table
	qry += "where s.icgc_mutation_id=g.icgc_mutation_id and s.pathogenic_estimate=1  "
	qry += "and s.reliability_estimate=1 "
	group_string = (",".join([quotify(gene2) for gene2 in gene_list]))
	qry += "and g.gene_symbol in (%s)" %  group_string

	ret = search_db(cursor,qry)
	if not ret:
		search_db(cursor,qry, verbose=True)
		exit()
	return ret[0][0]

########################################
def find_chromosome(cursor, gene):
	qry = "select chromosome from icgc.hgnc where approved_symbol = '%s'" % gene
	ret = search_db(cursor,qry)
	if not ret or ret ==[]:
		print("chromosome not found for %s (?)"%gene)
		search_db(cursor,qry,verbose=True)
		exit()
	return ret[0][0]

########################################
def get_donors(cursor, table):
	qry  = "select distinct(icgc_donor_id) from %s " % table
	return [ret[0] for ret in search_db(cursor,qry)]

def get_mutations(cursor, table):
	qry = "select  distinct(icgc_mutation_id)  from %s " % table
	return [ret[0] for ret in search_db(cursor,qry)]

def get_number_of_path_mutations_per_specimen(cursor, table, specimen_id):
	qry  = "select count( distinct icgc_mutation_id)  from %s " % table
	qry += "where icgc_specimen_id = '%s' " % specimen_id
	qry += "and pathogenic_estimate=1 and reliability_estimate=1"
	return search_db(cursor,qry)[0][0]

def get_consequence(cursor, chromosome, mutation):
	qry  = "select consequence, aa_mutation from mutations_chrom_%s " % chromosome
	qry += "where icgc_mutation_id='%s' " % mutation
	return search_db(cursor,qry)[0]

def get_specimens_from_donor(cursor, table, icgc_donor_id):
	qry = "select  distinct(icgc_specimen_id)  from %s " % table
	qry += "where icgc_donor_id = '%s'" % icgc_donor_id
	return [r[0] for r in search_db(cursor,qry)]

def get_specimen_type(cursor, tumor_short, spec_ids):
	specimen_type = {}
	for spec_id in spec_ids:
		qry = " select specimen_type from %s_specimen " % tumor_short
		qry += "where icgc_specimen_id = '%s'" % spec_id
		specimen_type[spec_id] = search_db(cursor,qry)[0][0]
	return specimen_type

def get_mutations_from_donor(cursor, table, icgc_donor_id):
	qry = "select  distinct(icgc_mutation_id)  from %s " % table
	qry += "where icgc_donor_id = '%s'" % icgc_donor_id
	return [r[0] for r in search_db(cursor,qry)]

def mutation_provenance(cursor, table, icgc_donor_id, icgc_mutation_id):
	qry = "select  distinct(icgc_specimen_id)  from %s " % table
	qry += "where icgc_donor_id='%s' and icgc_mutation_id='%s'" % (icgc_donor_id, icgc_mutation_id)
	return [r[0] for r in search_db(cursor,qry)]

#########################################
def mutations_in_gene_old(cursor, approved_symbol):
	qry  = "select ensembl_gene_id_by_hgnc, ensembl_gene_id, chromosome from hgnc "
	qry += "where approved_symbol = '%s'" % approved_symbol
	ensembl_gene_id_by_hgnc, ensembl_gene_id, chromosome = search_db(cursor,qry)[0]
	if ensembl_gene_id_by_hgnc != ensembl_gene_id:
		print("Ensembl id mismatch: (ensembl_gene_id_by_hgnc, ensembl_gene_id)")
		print(approved_symbol, ensembl_gene_id_by_hgnc, ensembl_gene_id)
		exit()

	qry  = "select m.icgc_mutation_id from mutations_chrom_%s m, locations_chrom_%s l "  % (chromosome, chromosome)
	qry += "where m.pathogenic_estimate=1 and m.start_position=l.position "
	qry += "and l.gene_relative like '%%%s%%' " % ensembl_gene_id
	ret = search_db(cursor,qry, verbose=True)
	if not ret: return []
	return [r[0] for r in ret]

#########################################
def mutations_in_gene(cursor, approved_symbol):
	qry  = "select icgc_mutation_id from mutation2gene "
	qry += "where gene_symbol='%s'" % approved_symbol
	ret = search_db(cursor,qry, verbose=True)
	if not ret: return []
	return [r[0] for r in ret]

#########################################
def pathogenic_mutations_in_gene(cursor, approved_symbol, chromosome, use_reliability=True):
	qry  = "select map.icgc_mutation_id from mutation2gene map, mutations_chrom_%s mut " % chromosome
	qry += "where map.gene_symbol='%s' " % approved_symbol
	qry += "and map.icgc_mutation_id=mut.icgc_mutation_id "
	qry += "and mut.pathogenic_estimate=1 "
	if use_reliability: qry += "and mut.reliability_estimate=1 "
	ret = search_db(cursor,qry, verbose=True)
	if not ret: return []
	return [r[0] for r in ret]

#########################################
def try_to_resolve(cursor, old_ensembl_gene_id):
	# not sure how stable or reliable this is, but if there is no common transcript,
	# something iss seriously foul with the old_ensembl_gene_id
	switch_to_db(cursor,"homo_sapiens_core_91_38")
	qry = "select distinct(gene_stable_id) from gene_archive "
	qry += "where transcript_stable_id in  "
	qry +="(select transcript_stable_id from gene_archive where gene_stable_id = '%s')" % old_ensembl_gene_id
	ret = search_db(cursor,qry)
	if not ret: return None

	candidate_ids = [r[0] for r in ret if r[0]!=old_ensembl_gene_id]
	if len(candidate_ids)==0: return None  # there should be another identifier, besides the one we started from

	latest_ensembl_entries = []
	for  candidate in candidate_ids:
		qry = "select * from gene  where stable_id = '%s'" % candidate
		if not search_db(cursor,qry): continue  # this is another old identifier
		latest_ensembl_entries.append(candidate)
	#I'm not sure what to make of this if there are two live identifiers
	# that the old one maps to
	if len(latest_ensembl_entries) != 1: return None

	new_ensembl_gene_id = latest_ensembl_entries[0]

	qry = "select approved_symbol from icgc.hgnc where ensembl_gene_id='%s'"% new_ensembl_gene_id
	ret = search_db(cursor,qry)
	if not ret: return None
	return ret[0][0]


def get_approved_symbol(cursor, ensembl_gene_id):
	qry = "select approved_symbol from hgnc where ensembl_gene_id='%s'"% ensembl_gene_id
	ret = search_db(cursor,qry)
	if not ret:
		symbol = try_to_resolve(cursor, ensembl_gene_id)
		# if it cannot be resolved, just use the ensembl_id_itself
		if not symbol: symbol=ensembl_gene_id
		switch_to_db(cursor,"icgc")
	else:
		symbol = ret[0][0]
	return symbol

#########################################
def gene_stable_id_2_canonical_transcript_id(cursor, gene_stable_id):
	qry  = "select  distinct(canonical_transcript) from ensembl_ids where  gene ='%s' " % gene_stable_id
	ret = search_db(cursor,qry)
	if not ret or len(ret) != 1:
		print("Warning: no unique canonical id could be found for %s" % gene_stable_id)
		return None
	return ret[0][0]


def list_of_transcript_ids_2_canonical_transcript_id(cursor, list_of_stable_transcript_ids):
	# list_od_stable_transcript_ids - refers to ensembl --> ENST00... identifier
	ensts = ",".join(["'%s'"%enst for enst in list_of_stable_transcript_ids])
	qry  = "select distinct(canonical_transcript) from ensembl_ids  where transcript in  (%s) " % ensts
	ret = search_db(cursor,qry)
	if not ret or len(ret) != 1:
		print("Warning: no unique canonical transcript id could be found for %s" % ensts)
		print("Qry was: ", qry)
		return None
	return ret[0][0]


