#! /usr/bin/python

# careful: mutation ids are unieq for position and mutation type, not for donor!

from icgc_utils.mysql   import  *
from icgc_utils.common_queries   import  *
import operator

verbose = False

#########################################
benign = ["downstream", "3_prime_UTR", "upstream", "intron", "synonymous", "5_prime_UTR",
	          "splice_region", "intragenic", "intergenic_region"]
disruptive = ["missense", "splice", "frame", "stop", "start"]

#
#=================
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

#=================
def mutation_description(cursor, table, icgc_donor_id, icgc_mutation_id):
	descr = {}
	symbol = {}
	qry  = "select consequence_type, gene_affected, transcript_affected from %s " % table
	qry += "where icgc_donor_id='%s' and icgc_mutation_id='%s'" % (icgc_donor_id, icgc_mutation_id)

	consequence_full = {}
	ret = search_db(cursor,qry)
	if not ret: return descr
	for consequence_type, gene_affected, transcript_affected in ret:
		if not gene_affected or gene_affected=="": continue
		if not consequence_full.has_key(gene_affected): consequence_full[gene_affected] = {}
		if not consequence_full[gene_affected] .has_key(transcript_affected): consequence_full[gene_affected][transcript_affected] = ""
		consequence_full[gene_affected][transcript_affected] += consequence_type
	for gene, cons_transcript in consequence_full.iteritems():
		symbol = get_approved_symbol(cursor, gene)
		for transcript, cons in cons_transcript.iteritems():
			dis = [d for d in disruptive if d in cons.lower()]
			if len(dis)==0: continue
			for d in dis:
				if not descr.has_key(symbol): descr[symbol]="";
				if not d in descr[symbol]: descr[symbol] += d+"; "

	return descr

#########################################
def mutation_characterization(cursor,  table, icgc_donor_id, mut_list, freq):
	genes_seen = set([])
	ct = 0
	for icgc_mutation_id in mut_list:
		ct +=1
		if ct%1000==0: print "    ***** ", ct
		descr = mutation_description(cursor, table, icgc_donor_id, icgc_mutation_id)
		if len(descr)==0: continue
		for gene, cons in descr.iteritems():
			if verbose: print "\t\t", icgc_mutation_id, gene, cons
			genes_seen.add(gene)
	# freq of mutations per patient
	for gene in genes_seen:
		if not freq.has_key(gene): freq[gene] = 0
		freq[gene] += 1

#########################################
#########################################
def main():


	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################

	tables = ['BRCA_simple_somatic']

	switch_to_db(cursor,"icgc")
	for table in tables:
		print "================================="
		print table
		# total number of donors?
		donors = get_donors(cursor, table)
		print "\t donors: ", len(donors)

		#==================
		#donors = ['DO891']
		freq_single_sample = {}
		freq_common_to_multiple_samples = {}

		ct = 0
		for icgc_donor_id in donors:
			ct += 1
			if ct%10==0: print " ** ", ct
			if verbose: print "\t", icgc_donor_id
			spec_ids  =  get_specimens_from_donor(cursor, table, icgc_donor_id)
			spec_type = get_specimen_type(cursor, table, spec_ids)
			if verbose:
				for spec_id, stype in spec_type.iteritems():
					print "\t\t", spec_id, stype
			mutation_ids = get_mutations_from_donor(cursor, table, icgc_donor_id)
			if verbose: print "\t\ttotal mutations: ",  len(mutation_ids)
			if len(spec_ids)==1:
				mutation_characterization(cursor,  table, icgc_donor_id, mutation_ids, freq_single_sample)
			else:
				muts_in_all_specs = []
				source = {}
				for spec_id in spec_ids: source[spec_id] = 0
				for icgc_mutation_id in mutation_ids:
					prov = mutation_provenance(cursor, table, icgc_donor_id, icgc_mutation_id)
					if len(prov) == len(spec_ids): muts_in_all_specs.append(icgc_mutation_id)
					for spec_id in prov:
						source[spec_id] += 1
				if verbose:
					print "\t\tin all specimens: ", len(muts_in_all_specs)
					for spec_id in spec_ids:
						print "\t\tin ", spec_id+"(%s):" % spec_type[spec_id][0], source[spec_id]
				mutation_characterization(cursor,  table, icgc_donor_id, muts_in_all_specs, freq_common_to_multiple_samples)
		freq = {k: freq_single_sample.get(k,0) + freq_common_to_multiple_samples.get(k,0)
		        for k in set(freq_single_sample) | set(freq_common_to_multiple_samples)}
		freq_sorted = sorted(freq.items(), key=operator.itemgetter(1))
		for gene, f in freq_sorted:
			if f<=1: continue
			print gene, f, freq_common_to_multiple_samples.get(gene,0)

		#exit()
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
