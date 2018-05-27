
from icgc_utils.mysql   import  *

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
def co_ocurrence_count(cursor, somatic_table, gene1, gene2):

	qry =  "select count(distinct s1.icgc_donor_id) ct "
	qry += "from mutation2gene g1, mutation2gene g2,  %s s1,  %s s2  " % (somatic_table, somatic_table)
	qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
	qry += "and s1.icgc_mutation_id=g1.icgc_mutation_id and g1.gene_symbol='%s' " % gene1
	qry += "and s2.icgc_mutation_id=g2.icgc_mutation_id and g2.gene_symbol='%s' " % gene2
	qry += "and s1.pathogenic_estimate=1 and s1.reliability_estimate=1  "
	qry += "and s2.pathogenic_estimate=1 and s2.reliability_estimate=1 "
	ret = search_db(cursor,qry)
	if not ret:
		search_db(cursor,qry, verbose=True)
		exit()
	return ret[0][0]

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

########################################
def find_chromosome(cursor, gene):
	qry = "select chromosome from hgnc where approved_symbol = '%s'" % gene
	ret = search_db(cursor,qry)
	if not ret or ret ==[]:
		print "chromosome not found for %s (?)"%gene
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

def get_specimen_type(cursor, table, spec_ids):
	specimen_type = {}
	for spec_id in spec_ids:
		qry = " select specimen_type from %s " % table.replace("simple_somatic","specimen")
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
		print "Ensembl id mismatch: (ensembl_gene_id_by_hgnc, ensembl_gene_id)"
		print approved_symbol, ensembl_gene_id_by_hgnc, ensembl_gene_id
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
def canonical_transcript_id_from_gene_stable_id(cursor, gene_stable_id):
	# do I have ensembl database locally?
	qry = "show databases like 'homo_sapiens_core%'"
	ret = search_db(cursor,qry)
	if not ret or not 'homo' in ret[0][0]:
		print "no database like homo_sapiens_core% available (to find canonical tr ids)"
		exit()
	ensembl_homo_sapiens_db = ret[0][0]
	qry  = "select t.stable_id from %s.gene g, %s.transcript t " % (ensembl_homo_sapiens_db, ensembl_homo_sapiens_db)
	qry += "where g.canonical_transcript_id=t.transcript_id "
	qry += "and g.stable_id='%s' " % gene_stable_id
	ret = search_db(cursor,qry)
	if not ret or len(ret) != 1:
		print "Warning: no unique canonical id could be found for %s" % gene_stable_id
		return None
	return ret[0][0]


def get_stable_id_for_canonical_transcript(cursor, list_of_stable_transcript_ids, gene_stable_id=None):

	# do I have ensembl database locally?
	qry = "show databases like 'homo_sapiens_core%'"
	ret = search_db(cursor,qry)
	if not ret or not 'homo' in ret[0][0]:
		print "no database like homo_sapiens_core% available (to resolve ENSTs)"
		exit()
	ensembl_homo_sapiens_db = ret[0][0]
	# list_od_stable_transcript_ids - refers to ensembl --> ENST00... identifier
	ensts = ",".join(["'%s'"%enst for enst in list_of_stable_transcript_ids])
	qry  = "select t.stable_id from %s.gene g, %s.transcript t " % (ensembl_homo_sapiens_db, ensembl_homo_sapiens_db)
	qry += "where g.gene_id = t.gene_id and  g.canonical_transcript_id=t.transcript_id "
	qry += "and t.stable_id in  (%s) " % ensts
	if gene_stable_id: # if we know which gene we are looking for, use that knowledge here
		# why all the fuss if I know the gene stable id?
		qry += "and g.stable_id='%s' " % gene_stable_id
	ret = search_db(cursor,qry)
	if not ret or len(ret) != 1:
		print "Warning: no unique canonical id could be found for %s" % ensts
		return None
	return ret[0][0]


