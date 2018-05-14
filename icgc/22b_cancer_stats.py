#! /usr/bin/python


from icgc_utils.common_queries   import  *


#########################################
def avg_number_of_muts_per_patient(cursor, table, donors, mutations_in_gene):
	number_of_patients_w_pathogenic_mutations = 0
	avg_no_muts = 0
	patients_with_muts_in_gene = {}
	for gene in mutations_in_gene.keys(): patients_with_muts_in_gene[gene] = 0
	if len(donors)==0: return
	for donor_id in  donors:
		qry  = "select distinct icgc_mutation_id from %s " % table
		qry += "where  pathogenic_estimate =1 "
		qry += "and   icgc_donor_id='%s' " % donor_id
		qry += "and  (total_read_count is null or (mutant_allele_read_count>3 and mut_to_total_read_count_ratio>0.2) )"
		ret = search_db(cursor,qry)
		number_of_muts = len(ret) if ret else 0
		if number_of_muts==0: continue

		number_of_patients_w_pathogenic_mutations += 1
		avg_no_muts += number_of_muts
		mut_ids = set([r[0] for r in ret])
		for gene, mutations in mutations_in_gene.iteritems():
			if len(mut_ids&set(mutations))>0: patients_with_muts_in_gene[gene] += 1

	if number_of_patients_w_pathogenic_mutations: avg_no_muts/= number_of_patients_w_pathogenic_mutations
	print "\t number of patients with pathogenic mutations: %d" % number_of_patients_w_pathogenic_mutations,
	print "\t (%d%%)" % ( float(number_of_patients_w_pathogenic_mutations)/len(donors)*100 )
	print "\t avg number of mutations  %.1f " %  avg_no_muts
	print "\t patients with mutations in "
	for gene, no_patients in  patients_with_muts_in_gene.iteritems():
		print "\t\t ", gene, no_patients


#########################################
#########################################

# produce table of the format
# tumor short | tumor long | number of patients | avg number of mutations per patient |
#  number of patients with mutated rpl5 (%of patients; number of genes which are seen mutated in the same or bigger number of patients)
#  | ditto for rp111

def main():

	db     = connect_to_mysql()
	cursor = db.cursor()
	#########################
	switch_to_db(cursor,"icgc")
	#########################

	gene_mutations = {}
	target_genes = ['RPL5']
	for gene in target_genes:
		chrom = find_chromosome(cursor, gene)
		print
		print "pathogenic_mutations in", gene, "chromosome", chrom
		gene_mutations[gene] = pathogenic_mutations_in_gene(cursor, gene, chrom)
		totmus = len(gene_mutations[gene])
		print "\t total mutations:", totmus
		print "\t the first 10:",  gene_mutations[gene][:10]
	print

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################



	for table in tables:
		#if table=="BRCA_simple_somatic": continue
		print "================================="
		print table

		# total number of donors?
		qry  = "select distinct(icgc_donor_id) from %s " % table
		donors = [ret[0] for ret in search_db(cursor,qry)]
		print "\t donors: ", len(donors)
		qry  = "select distinct(icgc_specimen_id) from %s " % table
		specimens = [ret[0] for ret in search_db(cursor,qry)]
		print "\t specimens: ", len(specimens)

		# number of unique mutations for each patient
		avg_number_of_muts_per_patient(cursor, table, donors, gene_mutations)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
'''
what are the missing mutations - the ones we dropped bcs of the low read quality?

	for mid in gene_mutations['RPL5']:
		for table in tables:
			qry = " select * from %s where icgc_mutation_id='%s'" %(table, mid)
			ret = search_db(cursor,qry)
			if not ret: continue
			print mid, table
			print qry
			#print "\n".join([",".join(r[0]) for r in ret])
			for r in ret:
				print "\t".join([str(field) for field in r])
			print
	exit()



'''