#! /usr/bin/python


import MySQLdb
from icgc_utils.mysql   import  *
#########################################
def get_canonical_transcript_ids(cursor, ensembl_gene_ids):
	switch_to_db(cursor, "homo_sapiens_core_91_38")
	qry = "select canonical_transcript_id from gene where stable_id in (%s)"%",".join(["'%s'"%g for g in ensembl_gene_ids])
	canonical_transcript_ids = [field[0] for field in  search_db(cursor,qry)]

	qry = "select stable_id from transcript where transcript_id in (%s)" % ",".join(["%d"%g for g in canonical_transcript_ids])
	stable_transcript_ids = [field[0] for field in  search_db(cursor,qry)]
	switch_to_db(cursor, "icgc")
	return  stable_transcript_ids


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

	#tables = ['EOPC_simple_somatic']

	db_name =  "icgc"
	switch_to_db(cursor, db_name)

	genes = ["TP53", "RPL5"]
	qry = "select ensembl_gene_id, approved_symbol from hgnc where approved_symbol in (%s)"%",".join(["'%s'"%g for g in genes])
	gene_name = {}

	for line in search_db(cursor,qry):
		[ensembl_gene_id, approved_symbol] = line
		gene_name[ensembl_gene_id] = approved_symbol
	ensembl_gene_ids = 	gene_name.keys()


	can_tr_ids = get_canonical_transcript_ids(cursor, ensembl_gene_ids)


	benign = ["downstream", "3_prime_UTR", "upstream", "intron", "synonymous", "5_prime_UTR",
	          "splice_region", "intragenic", "intergenic_region"]
	for table in tables:
		print "================================="
		print table
		# total number of donors?
		qry  = "select icgc_donor_id, submitted_sample_id, icgc_mutation_id, consequence_type, gene_affected, transcript_affected from %s " % table
		qry += "where gene_affected in (%s)  "%",".join(["'%s'"%g for g in ensembl_gene_ids])
		donors = {}
		donor_id_translation = {}
		for line in search_db(cursor,qry):
			[icgc_donor_id, submitted_sample_id, icgc_mutation_id, consequence_type, gene_affected, transcript_affected] = line
			donor_id_translation[icgc_donor_id] = submitted_sample_id
			if consequence_type==None or consequence_type=="" or consequence_type in benign: continue
			if not donors.has_key(icgc_donor_id): donors[icgc_donor_id] = {}
			if not donors[icgc_donor_id].has_key(gene_affected): donors[icgc_donor_id][gene_affected] = {}
			if not donors[icgc_donor_id][gene_affected].has_key(icgc_mutation_id):
				donors[icgc_donor_id][gene_affected][icgc_mutation_id] = {}
			if not donors[icgc_donor_id][gene_affected][icgc_mutation_id].has_key(transcript_affected):
				donors[icgc_donor_id][gene_affected][icgc_mutation_id][transcript_affected] = consequence_type
			else:
				donors[icgc_donor_id][gene_affected][icgc_mutation_id][transcript_affected] += "; " + consequence_type

		for donor, genes_affected in donors.iteritems():
			#print donor, genes_affected.keys()
			gn = [gene_name[e] for e in genes_affected.keys()]
			if 'RPL5' in gn:
				print donor, donor_id_translation[donor], [gene_name[e] for e in genes_affected.keys()]
			#for gene, mutations in genes_affected.iteritems():
				#print "\t", gene_name[gene]
				#for mutation, transcripts in mutations.iteritems():
				#	print "\t\t", mutation
				#	for transcr, consequence in transcripts.iteritems():
				#		print "\t\t\t", transcr, consequence




	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()

