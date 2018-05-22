#! /usr/bin/python


from icgc_utils.common_queries   import  *

verbose = False

#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./icgc_utils/kernprof.py -l 22a_cancer_stats.py
# python -m line_profiler 22a_cancer_stats.py.lprof
# see here https://github.com/rkern/line_profiler#line-profiler
# the reason I am using local kernprof.py is that I don't know where pip
# installed its version (if anywhere)
# @profile
#########################################
def avg_number_of_muts_per_patient(cursor, table, donors):
	number_of_patients_w_pathogenic_mutations = 0
	avg_no_muts = 0

	if len(donors)==0: return
	for donor_id in  donors:
		qry  = "select distinct icgc_mutation_id from %s " % table
		qry += "where  icgc_donor_id='%s' " % donor_id
		qry += "and  pathogenic_estimate =1 "
		qry += "and  (total_read_count is null or (mutant_allele_read_count>3 and mut_to_total_read_count_ratio>0.2) )"
		ret = search_db(cursor,qry)
		if not ret: continue
		number_of_muts = len(ret) if ret else 0
		if number_of_muts==0: continue

		number_of_patients_w_pathogenic_mutations += 1
		avg_no_muts += number_of_muts

	# why is this in this function?
	if number_of_patients_w_pathogenic_mutations: avg_no_muts/= number_of_patients_w_pathogenic_mutations
	return number_of_patients_w_pathogenic_mutations, avg_no_muts


#########################################
def patients_per_gene_breakdown(cursor, table):

	# this hinges on s.icgc_mutation_id=g.icgc_mutation_id
	# having icgc_mutation_id indexed both on s and g:
	qry  = "select g.gene_symbol symbol, count(distinct  s.icgc_donor_id) ct "
	qry += "from mutation2gene g, %s s  " % table
	qry += "where s.icgc_mutation_id=g.icgc_mutation_id and s.pathogenic_estimate=1  "
	qry += "and (s.total_read_count is null or (s.mutant_allele_read_count>3 and s.mut_to_total_read_count_ratio>0.2)) "
	qry += "group by symbol"
	ret = search_db(cursor,qry)
	if not ret:
		search_db(cursor,qry, verbose=True)
		exit()
	return dict(ret)

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
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")

	outf = open("mutations_per_cancer_breakdown.tsv", "w")
	outf.write("\t".join(["tumor short", "donors", "specimen",
							"donors w pathogenic mutations", "pct",
							"avg mutations per patient",
							"RPL5 donors", "pct of all donors",
							"genes wiht path muts in more donors than RPL5", "pct genome",
							"RPL11 donors", "pct of all donors",
							"genes wiht path muts in more donors than RPL11", "pct genome"
							])+"\n")

	#tables = ["BRCA_simple_somatic"]
	for table in tables:
		#if table=="BRCA_simple_somatic": continue

		tumor_short = table.split("_")[0]
		fields = [tumor_short]
		if verbose: print "================================="
		if verbose: print table

		# total number of donors?
		qry  = "select distinct(icgc_donor_id) from %s " % table
		donors = [ret[0] for ret in search_db(cursor,qry)]
		if verbose: print "\t donors: ", len(donors)
		fields.append(len(donors))

		qry  = "select distinct(icgc_specimen_id) from %s " % table
		specimens = [ret[0] for ret in search_db(cursor,qry)]
		if verbose: print "\t specimens: ", len(specimens)
		fields.append(len(specimens))

		# number of unique mutations for each patient
		number_of_patients_w_pathogenic_mutations,avg_no_muts = avg_number_of_muts_per_patient(cursor, table, donors)
		if verbose: print "\t number of patients with pathogenic mutations: %d" % number_of_patients_w_pathogenic_mutations,
		pct = float(number_of_patients_w_pathogenic_mutations)/len(donors)*100
		if verbose: print "\t (%d%%)" % (pct)
		fields.append(number_of_patients_w_pathogenic_mutations)
		fields.append("%.0f" % pct)

		if verbose: print "\t avg number of mutations  %.1f " % avg_no_muts
		fields.append(" %.1f " % avg_no_muts)

		patients_with_muts_in_gene = patients_per_gene_breakdown(cursor, table)

		if verbose: print "\t patients with mutations in "
		for gene in ['RPL5', 'RPL11']:
			if patients_with_muts_in_gene.has_key(gene):
				nr_muts =  patients_with_muts_in_gene[gene]
				genes_w_eq_or_gt_number_of_donors =  len([g for g in patients_with_muts_in_gene.keys() if patients_with_muts_in_gene[g]>=nr_muts])
				if verbose: print "\t\t ", gene, nr_muts,  nr_muts, genes_w_eq_or_gt_number_of_donors
				pct  = float(nr_muts)/number_of_patients_w_pathogenic_mutations*100
				fields.append(nr_muts)
				fields.append("%.1f" % pct)
				fields.append(genes_w_eq_or_gt_number_of_donors)
				fields.append("%.0f" % (genes_w_eq_or_gt_number_of_donors/20000.0*100))
			else:
				if verbose: print "\t\t ", gene,0
				fields.append(0)
				fields.append(0.0)
				fields.append(20000)
				fields.append(100)
		outf.write("\t".join([str(f) for f in fields])+"\n")
		outf.flush()

	outf.close()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
