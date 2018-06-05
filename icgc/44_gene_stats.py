#! /usr/bin/python


from icgc_utils.common_queries   import  *

verbose =True

#########################################
# sample codes:
# 01	Primary solid Tumor	TP
# 02	Recurrent Solid Tumor	TR
# 03	Primary Blood Derived Cancer - Peripheral Blood	TB
# 04	Recurrent Blood Derived Cancer - Bone Marrow	TRBM
# 05	Additional - New Primary	TAP
# 06	Metastatic	TM
# 07	Additional Metastatic	TAM
# 08	Human Tumor Original Cells	THOC
# 09	Primary Blood Derived Cancer - Bone Marrow	TBM
# 10	Blood Derived Normal	NB
# 11	Solid Tissue Normal	NT
# 12	Buccal Cell Normal	NBC
# 13	EBV Immortalized Normal	NEBV
# 14	Bone Marrow Normal	NBM
# 20	Control Analyte	CELLC
# 40	Recurrent Blood Derived Cancer - Peripheral Blood	TRB
# 50	Cell Lines	CELL
# 60	Primary Xenograft Tissue	XP
# 61	Cell Line Derived Xenograft Tissue	XCL


def spec_from_TCGA(sample_barcode):
	# we want to translate this to something similar to what ICGC is using
	# roughly: Normal, Primary, Metastatic, Recurrent, Cell_line
	tcga_sample_code = sample_barcode.split("-")[3][:2]
	if tcga_sample_code in ['01','03','05', '09']:
		return 'Primary'

	elif tcga_sample_code in ['02','04','40']:
		return 'Recurrent'

	elif tcga_sample_code in ['06','07']:
		return 'Metastatic'

	elif tcga_sample_code in ['10','11','12','13','14']:
		return 'Normal'

	elif tcga_sample_code in ['08','50']:
		return 'Cell_line'

	return "Other"




#########################################
def gene_mutations(cursor, tumor_short, gene):

	qry  = "select m.icgc_donor_id, m.submitted_sample_id, m.chromosome, m.control_genotype, m.tumor_genotype, "
	qry += "g.icgc_mutation_id, m.icgc_specimen_id "
	qry += "from mutation2gene g,  %s_simple_somatic m " % tumor_short
	qry += "where g.gene_symbol='%s' " % gene
	qry += "and g.icgc_mutation_id = m.icgc_mutation_id "
	qry += "and m.pathogenic_estimate=1 and m.reliability_estimate=1"

	ret = search_db(cursor,qry)
	if not ret: return

	donor_rows = {}
	p53_status_per_specimen = {}
	specimen_seen = {}
	mutations_per_specimen = {}
	for line in ret:
		[donor, sample, chromosome,  cgenotype, tgenotype, mutation, specimen] = line

		if donor[2]=="T":
			# this came from TCGA => we do not have the specimen table,
			# but we can tell its type from the TCGA code itself
			spec_type = spec_from_TCGA(sample)
		else:
			# find the specimen type from the specimen table
			spec_type = get_specimen_type(cursor,"%s_specimen"%tumor_short,[specimen])[specimen]

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

		if not specimen: specimen=""
		if not cgenotype: cgenotype=""
		if not donor_rows.has_key(donor):
			entry = "\t".join([tumor_short, donor, specimen, spec_type[:1], no_muts, mutation,  cgenotype, tgenotype, consequence, aa_change, out_p53_status])
			donor_rows[donor] = [entry]
		else:
			entry = "\t".join([tumor_short, "", specimen, spec_type[:1], no_muts, mutation,  cgenotype, tgenotype,  consequence, aa_change, out_p53_status])
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
