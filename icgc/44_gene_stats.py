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
def gnomad_info(cursor, mutation, chromosome):

	# gnomad has info for standard chromosomes and X
	if chromosome not in [str(i) for i in range(23)] + ['X']: return "n/a"

	returnval =  "n/a"
	# check genome assembly:
	# one more time:
	# GRCh38 = hg38
	# GRCh37 = hg19 <--- gnomad
	# NCBI Build 36.1 = hg18
	[start_position, assembly, reference_genome_allele, mutated_to_allele] = [None]*4
	qry  = "select start_position, assembly, reference_genome_allele, mutated_to_allele "
	qry += "from mutations_chrom_%s " % chromosome
	qry += "where icgc_mutation_id='%s'" % mutation

	ret = search_db(cursor,qry,verbose=True)
	if len(ret)>1:
		print "we should not be here: multiple returns for \n", qry
		exit()
	if ret:
		[start_position, assembly, reference_genome_allele, mutated_to_allele] = ret[0]
		if not assembly in ['GRCh37','hg19']:
			# gnomad chromosome loataions use this assembly
			# TODO: one day - translation to other assemblies (?)
			assembly = None
	else:
		assembly = None

	if assembly:
		qry  = "select variant_count,total_count  from gnomad.gnomad_freqs_chr_%s " % chromosome
		qry += "where position=%d " % start_position
		qry += "and reference='%s' and variant='%s' " % (reference_genome_allele,mutated_to_allele )

		ret = search_db(cursor,qry)
		if ret:
			if len(ret)>1:
				print "we should not be here: multiple returns for \n", qry
				exit()
			if type(ret[0][0]) in [int,long] and type(ret[0][1]) in [int,long]:
				[variant_count,total_count] = ret[0]
				returnval = "%.1e"%(float(variant_count)/total_count)
			else:
				print "unexpected type for variant and total counts"
				print search_db(cursor,qry,verbose=True)
				exit()
		else:
			returnval = "0.0"
	return returnval

patient_count = 0

#########################################
def donor_mutations_to_printable_format(cursor, tumor_short, donor_mutations, exp_results, hide_id=True):

	p53_status_per_specimen = {}
	mutations_per_specimen  = {}
	specimen_seen = {}
	donor_rows = []
	donor_ct = 0
	for donor, mutations in donor_mutations.iteritems():
		global patient_count
		patient_count += 1
		donor_display_name = donor
		for mutation, mtn_info in mutations.iteritems():

			chromosome = mtn_info['chromosome']
			freq_in_gen_population = gnomad_info(cursor, mutation, chromosome)

			###################
			# sample types (primary, metastatic etc)
			specimen_type_short = []
			for sample in mtn_info['samples']:
				idx = mtn_info['samples'].index(sample)
				specimen = mtn_info['specimens'][idx]
				if sample[:4] in ['TCGA','TARG']:
					# this came from TCGA => sometimes we do not have the specimen table,
					# but we can tell its type from the TCGA code itself
					spec_type = spec_from_TCGA(sample)
				else:
					# find the specimen type from the specimen table
					spec_type = get_specimen_type(cursor,tumor_short,[specimen])[specimen]
				specimen_type_short.append(spec_type[0])

				if sample[:4]=='TCGA': donor_display_name = sample[0:15] # for comparison with previous results

			###################
			# specimen related info
			p53_gist   = "wt"
			p53_detail = ""
			specimen_number_of_mutations = []
			for specimen in mtn_info['specimens']:
				if not specimen_seen.has_key(specimen): # we might have seen it related to another mutation
					specimen_seen[specimen] = True
					table = "%s_simple_somatic" % tumor_short
					mutations_per_specimen[specimen]  = get_number_of_path_mutations_per_specimen(cursor, table, specimen)
					p53_status_per_specimen[specimen] = find_53_status(cursor, tumor_short, specimen)

				specimen_number_of_mutations.append(str(mutations_per_specimen[specimen]))
				if p53_gist!="pathogenic":
					new_gist, detail = p53_status_per_specimen[specimen]
					if new_gist!="wt":
						p53_gist=new_gist
						p53_detail = detail

			###################
			# output row for this mutation

			consequence, aa_change = get_consequence(cursor, chromosome, mutation)
			if not consequence: consequence = ""
			aa_change = aa_change_cleanup(cursor, aa_change)
			exp = exp_results.get(aa_change,"")

			if consequence=='frameshift':
				exp = " | ".join(['x']*5)

			if exp == " | ".join(['o']*5): exp = ""

			if hide_id: donor_display_name=str(patient_count)
			entry = "\t".join([tumor_short, donor_display_name,  ",".join(specimen_type_short), ",".join(specimen_number_of_mutations),
			                   mtn_info['cgenotype'], mtn_info['tgenotype'], consequence,
				               aa_change, freq_in_gen_population, p53_gist, p53_detail , exp])
			entry = entry.replace("_"," ")
			donor_rows.append(entry)

	return "\n".join(donor_rows)


#########################################
def gene_mutations(cursor, table, gene):

	qry  = "select m.icgc_donor_id, m.submitted_sample_id, m.chromosome, m.control_genotype, m.tumor_genotype, "
	qry += "g.icgc_mutation_id, m.icgc_specimen_id "
	qry += "from mutation2gene g,  %s m " % table
	qry += "where g.gene_symbol='%s' " % gene
	qry += "and g.icgc_mutation_id = m.icgc_mutation_id "
	qry += "and m.pathogenic_estimate=1 and m.reliability_estimate=1"

	ret = search_db(cursor,qry)
	if not ret: return

	donor_mutations = {}

	for line in ret:
		[donor, sample, chromosome,  cgenotype, tgenotype, mutation, specimen] = line
		if not donor_mutations.has_key(donor):
			donor_mutations[donor] = {}
		if not donor_mutations[donor].has_key(mutation):
			donor_mutations[donor][mutation] = {
				'chromosome' : chromosome,
				'cgenotype' : cgenotype if cgenotype else "",
				'tgenotype' : tgenotype if tgenotype else "",
				'specimens' : [specimen],
				'samples'   : [sample]
			}
		else:
			if donor_mutations[donor][mutation]['chromosome'] != chromosome:
				print "different chromosome for mutation %s (?!)" % mutation
				exit(1)
			if donor_mutations[donor][mutation]['cgenotype'] !="" and cgenotype !="" and \
					donor_mutations[donor][mutation]['cgenotype'] != cgenotype:
				print "different cgenotype for mutation %s " % mutation
				print "in ", donor_mutations[donor][mutation]['specimens'], "and", specimen
				print donor_mutations[donor][mutation]['cgenotype'], "vs", cgenotype
				exit(1)
			if donor_mutations[donor][mutation]['tgenotype'] !="" and tgenotype !="" and \
					donor_mutations[donor][mutation]['tgenotype'] != tgenotype:
				print "different tgenotype for mutation %s " % mutation
				print "in ", donor_mutations[donor][mutation]['specimens'], "and", specimen
				print donor_mutations[donor][mutation]['tgenotype'], "vs", tgenotype
				exit(1)

			donor_mutations[donor][mutation]['specimens'].append(specimen)
			donor_mutations[donor][mutation]['samples'].append(sample)


	return donor_mutations

#########################################
def parse_exp(exp_results_file):
	# I am assuming that the columns in this file look like this:
	# MUTATION	; Tumor type;	p53 status;	5S rRNA binding;	Mdm2/RPL5/RPL11complex formation; \
	# aggregation;	p53 activation;	incorporation into ribosomes
	translation = {"YES":"+", "NO":"-", "NT": "o", "P":"p"}
	exp_results_hash = {}
	inf = open(exp_results_file,"r")
	for line in inf:
		field =[f.replace(" ","") for f in  line.rstrip().split("\t")]
		if 'mutation' in field[0].lower():
			continue
		else:
			# turn aggregation into 'proper folding'
			# or some such, so that "+" means "function unaffected"
			if field[5] == "YES":
				field[5]="NO"
			elif field[5] == "NO":
				field[5]="YES"
			if exp_results_hash.has_key(field[0]):
				if exp_results_hash[field[0]] != " | ".join([translation[val] for val in field[3:]]):
					print "Note: result mismatch for", field[0], " in ", exp_results_file
					#print exit()
			# extra: if the mutant does not fold, it presumably has no other function
			if field[5]=="NO":
				exp_results_hash[field[0]] = " | ".join(['x','x','-','x','x'])
			else:
				exp_results_hash[field[0]] = " | ".join([translation[val] for val in field[3:]])
	inf.close()

	return exp_results_hash

#########################################
#########################################
# produce table of the format
# donor id  #muts_in_sample  control_genotype  tumor_genotype  aa_change1     gene2     variant2   homo/hetero    aa_change2

def main():

	gene = 'RPL11'

	exp_results_file = "/home/ivana/Dropbox/Sinisa/ribosomal/rezultati/Ines_rezultati_Feb2018/rpl5.csv"
	if not os.path.exists(exp_results_file):
		print exp_results_file, "not found"
		exit()

	exp_results = parse_exp(exp_results_file)

	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")

	outf = open("%s_per_cancer_breakdown.tsv"%gene, "w")
	outf.write("\t".join(["tumor short",  "donor", "spec type",
							"no muts",   "control genotype", "tumor genotype",
							"consequence", "aa change", "freq in general population",
                            "p53 status", "p53 mutation","function"])+"\n")

	for table in tables:
		tumor_short = table.split("_")[0]
		if verbose: print "================================="
		if verbose: print table
		donor_mutations = gene_mutations(cursor, table, gene)
		if not donor_mutations: continue
		donor_rows = donor_mutations_to_printable_format(cursor, tumor_short, donor_mutations, exp_results, hide_id=True);
		if not donor_rows: continue
		outf.write(donor_rows+"\n")

	outf.close()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
