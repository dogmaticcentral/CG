#! /usr/bin/python


from icgc_utils.common_queries   import  *

verbose = False

#########################################
def find_53_status(cursor, tumor_short, specimen):
	qry  = "select g.icgc_mutation_id, m.pathogenic_estimate "
	qry += "from mutation2gene g, %s_simple_somatic m " % (tumor_short)
	qry += "where g.gene_symbol='TP53' "
	qry += "and m.icgc_specimen_id = '%s' "  % specimen
	qry += "and g.icgc_mutation_id = m.icgc_mutation_id"
	ret = search_db(cursor,qry)
	if not ret: return "wt"
	for line in ret:
		if line[1]==1: return "pathogenic"
	else:
		return "benign"


#########################################
def gene_mutations(cursor, tumor_short, gene):
	qry  = "select m.icgc_donor_id, m.submitted_sample_id, g.icgc_mutation_id, m.icgc_specimen_id, s.specimen_type "
	qry += "from mutation2gene g,  %s_simple_somatic m, %s_specimen s " % (tumor_short, tumor_short)
	qry += "where g.gene_symbol='%s' " % gene
	qry += "and g.icgc_mutation_id = m.icgc_mutation_id "
	qry += "and m.icgc_specimen_id = s.icgc_specimen_id "
	qry += "and m.pathogenic_estimate=1"

	ret = search_db(cursor,qry)
	if not ret: return
	muts_per_donor = {}
	p53_status_per_specimen = {}
	for line in ret:
		[donor, sample, mutation, specimen, spec_type] = line
		if sample[:4]=='TCGA': donor = sample[0:15]
		if not p53_status_per_specimen.has_key(specimen):
			p53_status_per_specimen[specimen] = find_53_status(cursor, tumor_short, specimen)
			out_p53_status = p53_status_per_specimen[specimen]
		else:
			out_p53_status = ""
		if not muts_per_donor.has_key(donor):
			entry = "\t".join([tumor_short, donor, mutation, specimen, spec_type[:1], out_p53_status])
			muts_per_donor[donor] = [entry]
		else:
			entry = "\t".join([tumor_short, "", mutation, specimen, spec_type[:1], out_p53_status])
			muts_per_donor[donor].append(entry)

	for donor, entries in muts_per_donor.iteritems():
		for entry in entries:
			print entry

#########################################
#########################################
# produce table of the format
# donor id  #muts_in_sample    gene1   variant1    homo/hetero  change1     gene2     variant2   homo/hetero    aa_change2

def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")

	outf= open ("RPL5_per_cancer_breakdown.tsv", "w")
	outf.write( "\t".join(["tumor short",])+"\n")

	#tables = ["BLCA_simple_somatic"]
	for table in tables:
		tumor_short = table.split("_")[0]
		fields = [tumor_short]
		if verbose: print "================================="
		if verbose: print table
		gene_mutations(cursor, tumor_short, 'RPL11')

	outf.close()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
