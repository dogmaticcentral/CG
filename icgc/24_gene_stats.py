#! /usr/bin/python


from icgc_utils.common_queries   import  *

verbose =True

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
def transcript_location_cleanup(cursor, loc, gene_stable_id):
	if not loc: return ""
	if loc== "": return loc
	location = {};
	for enst_loc in loc.split(";"):
		[e, c] = enst_loc.split(":")
		location[e] = c
	enst_canonical = get_stable_id_for_canonical_transcript(cursor, location.keys(), gene_stable_id)
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
	enst_canonical = get_stable_id_for_canonical_transcript(cursor, change.keys())
	if not enst_canonical: return aa_change
	return change[enst_canonical]


#########################################
def gene_mutations(cursor, tumor_short, gene):

	qry  = "select m.icgc_donor_id, m.submitted_sample_id, m.chromosome, m.control_genotype, m.tumor_genotype, "
	qry += "g.icgc_mutation_id, m.icgc_specimen_id, s.specimen_type "
	qry += "from mutation2gene g,  %s_simple_somatic m, %s_specimen s " % (tumor_short, tumor_short)
	qry += "where g.gene_symbol='%s' " % gene
	qry += "and g.icgc_mutation_id = m.icgc_mutation_id "
	qry += "and m.icgc_specimen_id = s.icgc_specimen_id "
	qry += "and m.pathogenic_estimate=1 and m.reliability_estimate=1"

	ret = search_db(cursor,qry)
	if not ret: return
	donor_rows = {}
	p53_status_per_specimen = {}
	specimen_seen = {}
	mutations_per_specimen = {}
	for line in ret:
		[donor, sample, chromosome,  cgenotype, tgenotype, mutation, specimen, spec_type] = line
		if sample[:4]=='TCGA': donor = sample[0:15]

		###################
		# specimen related info
		if not specimen_seen.has_key(specimen):
			specimen_seen[specimen] = True
			p53_status_per_specimen[specimen] = find_53_status(cursor, tumor_short, specimen)
			table = "%s_simple_somatic" % tumor_short
			mutations_per_specimen[specimen] = get_number_of_path_mutations_per_specimen(cursor, table, specimen)

		out_p53_status = "\t".join(p53_status_per_specimen[specimen])
		no_muts        = str(mutations_per_specimen[specimen])
		consequence, aa_change = get_consequence(cursor, chromosome, mutation)
		if not consequence: consequence = ""
		aa_change = aa_change_cleanup(cursor, aa_change)

		if not donor_rows.has_key(donor):
			entry = "\t".join([tumor_short, donor,specimen, spec_type[:1], no_muts, mutation,  cgenotype, tgenotype , consequence, aa_change, out_p53_status])
			donor_rows[donor] = [entry]
		else:
			entry = "\t".join([tumor_short, "",specimen, spec_type[:1], no_muts, mutation,  cgenotype, tgenotype,  consequence, aa_change, out_p53_status])
			donor_rows[donor].append(entry)

	return donor_rows



#########################################
#########################################
# produce table of the format
# donor id  #muts_in_sample  control_genotype  tumor_genotype  aa_change1     gene2     variant2   homo/hetero    aa_change2

def main():

	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")

	gene = 'RPL5'
	outf = open("%s_per_cancer_breakdown.tsv"%gene, "w")
	outf.write("\t".join(["tumor short",  "donor","specimen", "spec_type",
							"no_muts", "mutation",  "cancer_genotype", "tumor_genotype",
							"consequence", "aa_change", "p53_status"])+"\n")

	for table in tables:
		tumor_short = table.split("_")[0]
		fields = [tumor_short]
		if verbose: print "================================="
		if verbose: print table
		donor_rows = gene_mutations(cursor, tumor_short, gene)
		if not donor_rows: continue
		for donor, entries in donor_rows.iteritems():
			for entry in entries:
				outf.write(entry+"\n")

	outf.close()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
