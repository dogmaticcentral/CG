#! /usr/bin/python
import subprocess
import time

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *

tcga_icgc_table_correspondence = {
"ACC_somatic_mutations" : None,
"ALL_somatic_mutations" : "ALL_simple_somatic",
"BLCA_somatic_mutations": "BLCA_simple_somatic",
"BRCA_somatic_mutations": "BRCA_simple_somatic",
"CESC_somatic_mutations": "CESC_simple_somatic",
"CHOL_somatic_mutations": None,
"COAD_somatic_mutations": "COCA_simple_somatic",
"DLBC_somatic_mutations": "DLBC_simple_somatic",
"ESCA_somatic_mutations": "ESAD_simple_somatic",
"GBM_somatic_mutations" : "GBM_simple_somatic",
"HNSC_somatic_mutations": "HNSC_simple_somatic",
"KICH_somatic_mutations": "KICH_simple_somatic",
"KIRC_somatic_mutations": "KIRC_simple_somatic",
"KIRP_somatic_mutations": "KIRP_simple_somatic",
"LAML_somatic_mutations": "AML_simple_somatic",
"LGG_somatic_mutations" : "LGG_simple_somatic",
"LIHC_somatic_mutations": "LICA_simple_somatic",
"LUAD_somatic_mutations": "LUAD_simple_somatic",
"LUSC_somatic_mutations": "LUSC_simple_somatic",
"MESO_somatic_mutations": None,
"OV_somatic_mutations"  : "OV_simple_somatic",
"PAAD_somatic_mutations": "PACA_simple_somatic",
"PCPG_somatic_mutations": None,
"PRAD_somatic_mutations": "PRAD_simple_somatic",
"READ_somatic_mutations": "COCA_simple_somatic",
"SARC_somatic_mutations": "SARC_simple_somatic",
"SKCM_somatic_mutations": "MELA_simple_somatic",
"STAD_somatic_mutations": "GACA_simple_somatic",
"TGCT_somatic_mutations": None,
"THCA_somatic_mutations": "THCA_simple_somatic",
"THYM_somatic_mutations": None,
"UCEC_somatic_mutations": "UCEC_simple_somatic",
"UCS_somatic_mutations" : "UTCA_simple_somatic",
"UVM_somatic_mutations" : None
}

variant_columns = ['icgc_mutation_id', 'chromosome','icgc_donor_id', 'icgc_specimen_id', 'icgc_sample_id',
                   'submitted_sample_id','control_genotype', 'tumor_genotype', 'total_read_count', 'mutant_allele_read_count']

# we'll take care of 'aa_mutation' and 'consequence_type will be handled separately
mutation_columns = ['icgc_mutation_id', 'start_position', 'end_position', 'mutation_type',
					'mutated_from_allele', 'mutated_to_allele', 'reference_genome_allele']


location_columns = ['position', 'gene_relative', 'transcript_relative']

################################################################
# stop_retained: A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
consequence_vocab = ['stop_lost', 'synonymous', 'inframe_deletion', 'inframe_insertion', 'stop_gained',
                     '5_prime_UTR_premature_start_codon_gain',
                     'start_lost', 'frameshift', 'disruptive_inframe_deletion', 'stop_retained',
                     'exon_loss', 'disruptive_inframe_insertion', 'missense']

# location_vocab[1:4] is gene-relative
# location_vocab[1:4] is transcript-relative
location_vocab = ['intergenic_region', 'intragenic', 'upstream', 'downstream',
                  '5_prime_UTR', 'exon',  'coding_sequence', 'initiator_codon',
                  'splice_acceptor', 'splice_region', 'splice_donor',
                  'intron', '3_prime_UTR', ]

# this is set literal
pathogenic = {'stop_lost', 'inframe_deletion', 'stop_gained', '5_prime_UTR_premature_start_codon_gain',
                     'start_lost', 'frameshift', 'disruptive_inframe_deletion',
                     'exon_loss', 'disruptive_inframe_insertion', 'missense',
                     'splice_acceptor', 'splice_region', 'splice_donor'
             }


def create_icgc_table(cursor, tcga_table):

	icgc_table = tcga_table.split("_")[0]+"_simple_somatic"
	if check_table_exists(cursor, 'icgc', icgc_table): return

	switch_to_db(cursor,'icgc')

	qry  = ""
	qry += "  CREATE TABLE  %s (" % icgc_table
	qry += "     id INT NOT NULL AUTO_INCREMENT, "
	qry += "  	 icgc_mutation_id VARCHAR (20) NOT NULL, "
	qry += "     chromosome CHAR(2) NOT NULL,"
	qry += "  	 icgc_donor_id VARCHAR (20) NOT NULL, "
	qry += "     icgc_specimen_id VARCHAR (20), "
	qry += "     icgc_sample_id VARCHAR (20), "
	qry += "     submitted_sample_id VARCHAR (50), "
	qry += "	 control_genotype VARCHAR (430) NOT NULL, "
	qry += "	 tumor_genotype VARCHAR (430) NOT NULL, "
	qry += "     total_read_count INT, "
	qry += "     mutant_allele_read_count INT, "
	qry += "     mut_to_total_read_count_ratio float default 0.0,"
	qry += "     pathogenic_estimate boolean default 0,"

	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print qry
	print rows


#########################################
def quotify(something):
	if not something:
		return ""
	if type(something)==str:
		return "\'"+something+"\'"
	else:
		return str(something)




#########################################
def insert (cursor, table, columns, values):

	nonempty_values = []
	corresponding_columns = []
	for i in range(len(values)):
		if not values[i] or  values[i] == "": continue
		nonempty_values.append(values[i])
		corresponding_columns.append(columns[i])
	qry = "insert into %s (%s) " %(table, ",".join(corresponding_columns))
	qry += "values (%s) " % ",".join(nonempty_values)
	search_db(cursor, qry)



#########################################
def output_annovar_input_file (cursor, db_name, tcga_table, already_deposited_samples):

	switch_to_db(cursor,db_name)
	meta_table_name = tcga_table.split("_")[0] + "_mutations_meta"

	# which assemblies do I have here?
	# I'm counting on the assemblies to be hg18 and hg19 - I have annovar tables for those
	# otherwise the annovar should complain
	qry  = "select distinct m.assembly "
	qry += "from %s s, %s m " % (tcga_table, meta_table_name)
	qry += "where s.meta_info_id=m.id"
	assemblies = [ret[0] for ret in search_db(cursor, qry)]

	# get the info that annovar needs
	qry  = "select s.tumor_sample_barcode, s.chromosome, s.start_position, s.end_position, "
	qry += "s.reference_allele, s.tumor_seq_allele1, s.tumor_seq_allele2, m.assembly  "
	qry += "from %s s, %s m " % (tcga_table, meta_table_name)
	qry += "where s.meta_info_id=m.id"

	rows = search_db(cursor, qry)
	outf = {}
	outfname = {}
	for assembly in assemblies:
		outfname[assembly] = "%s.%s.avinput" % (tcga_table,assembly)
		outf[assembly] = open(outfname[assembly], 'w')

	for row in rows:
		(tumor_sample_barcode, chromosome, start_position, end_position, reference_allele,
			tumor_seq_allele1, tumor_seq_allele2, assembly) = row
		if tumor_sample_barcode in already_deposited_samples: continue
		differing_allele = tumor_seq_allele1
		if differing_allele==reference_allele: differing_allele = tumor_seq_allele2
		# somebody in  the TCGA decided toinnovat and mark the insert with the before- and after- postion
		# Annovar expects both numbers to be the same in such case
		if reference_allele=="-": end_position=start_position
		outrow = "%s\t%d\t%d\t%s\t%s" % (chromosome, start_position, end_position, reference_allele, differing_allele)
		print >> outf[assembly], outrow

	for assembly in assemblies:
		outf[assembly].close()
		subprocess.call(["bash","-c", "sort %s | uniq > %s.tmp" % (outfname[assembly], outfname[assembly])])
		os.rename(outfname[assembly]+".tmp", outfname[assembly])

	return outfname

##################################
def run_annovar(avinput, table_name):
	avoutname = {}
	for assembly, avinfile in avinput.iteritems():
		avoutname[assembly] = "%s.%s_multianno.txt" % (table_name, assembly)
		if not os.path.exists(avoutname[assembly]) or os.path.getsize(avoutname[assembly])==0:
			cmd  = "/home/ivana/third/annovar/table_annovar.pl %s " % avinfile
			cmd += "/home/ivana/third/annovar/humandb/ -buildver %s -out %s " % (assembly, table_name)
			cmd += " -protocol ensGene  -operation g  -nastring ."
			subprocess.call(cmd, shell=True)
			# clean the junk
			cmd = "rm %s.ensGene.variant_function " % table_name
			cmd +="%s.ensGene.exonic_variant_function %s.ensGene.log" % (table_name, table_name)
			subprocess.call ( cmd, shell=True)
	return avoutname

#########################################
def store_annotation (cursor, icgc_table, avoutput):
	# store location and mutation info
	# the do another sweep through the tcga table to associate with variants
	# TODO pwd
	# HERE !!!
	for assembly, avfile in avoutput.iteritems():
		inf = open (avfile, "r")
		for line in inf:
			if line[:3]=="Chr": continue
			fields = line.rstrip().split('\t')
			[chrom, start, end] = fields[:3]
			fields = fields[-1].split(',')[0].split(':')

			if len(fields)<2: continue
			[cdna_change_position, val1, val2] = parse_mutation(fields[-2])
			if not cdna_change_position: continue
			[aa_change_position, val1, val2] = parse_mutation(fields[-1].replace('p.','').replace(' ', ''))
			if not aa_change_position: continue
			aa_change = val1 + str(aa_change_position) + val2
			if val1==val2:
				classf = "silent"
			else:
				classf = "missense_mutation"


#########################################
def process_table(home, cursor, tcga_table, icgc_table, already_deposited_samples):

	# make a workdir and move there
	tumor_short = tcga_table.split("_")[0]
	workdir  = tumor_short + "_annovar"
	workpath = "{}/annovar/{}".format(home,workdir)
	if not os.path.exists(workpath): os.makedirs(workpath)
	os.chdir(workpath)

	# use tcga entries to create input for annovar
	avinput  = output_annovar_input_file (cursor, 'tcga', tcga_table, already_deposited_samples)
	for assembly, avinfile  in avinput.iteritems():
		print "{}/{}".format(workpath,avinfile)

	# fork annovar process
	avoutput = run_annovar(avinput,tcga_table)

	# store the annotated input to icgc tables - *_simple_somatic, mutations_chrom*, and locations_chrom_*
	# store_annotation (cursor, icgc_table, avoutput)


#########################################
def add_tcga_diff(tcga_tables, other_args):

	db     = connect_to_mysql()
	cursor = db.cursor()
	home = os.getcwd()
	for tcga_table in tcga_tables:

		icgc_table =  tcga_icgc_table_correspondence[tcga_table]

		# where in the icgc classification does this symbol belong?
		time0 = time.time()
		print
		print "===================="
		print "processing tcga table ", tcga_table, os.getpid()
		print "will be stored in ", icgc_table


		if not icgc_table: create_icgc_table(cursor,tcga_table)

		#tcga samples in tcga_table
		qry = "select distinct(tumor_sample_barcode) from tcga.%s" % tcga_table
		tcga_tumor_sample_ids = [ret[0] for ret in search_db(cursor,qry)]

		#tcga samples already deposited in icgc
		qry  = "select distinct(submitted_sample_id) from icgc.%s " % icgc_table
		qry += "where submitted_sample_id like 'tcga%'"
		ret = search_db(cursor,qry)
		icgc_tumor_sample_ids = [r[0] for r in ret] if ret else []

		#tcga samples already deposited in icgc
		already_deposited_samples = list(set(icgc_tumor_sample_ids).difference(set(tcga_tumor_sample_ids)))
		samples_not_deposited     = list(set(tcga_tumor_sample_ids).difference(set(icgc_tumor_sample_ids)))

		# I am taking a leap of faith here, and I believe that the deposited data is
		# really identical to what we have in tcga
		print "not deposited:", len(samples_not_deposited)

		process_table(home, cursor, tcga_table, icgc_table, already_deposited_samples)
		print "\t overall time for %s: %.3f mins; pid: %d" % (tcga_table, float(time.time()-time0)/60, os.getpid())


	cursor.close()
	db.close()

	return


#########################################
#########################################
def main():

	# divide by cancer types, because I have duplicates within each cancer type
	# that I'll resolve as I go, but I do not want the threads competing)
	db     = connect_to_mysql()
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='tcga' and table_name like '%_somatic_mutations'"
	tcga_tables = [field[0] for field in search_db(cursor,qry)]

	number_of_chunks = 1  # myISAM does not deadlock
	parallelize(number_of_chunks, add_tcga_diff, tcga_tables, [])

#########################################
if __name__ == '__main__':
	main()
