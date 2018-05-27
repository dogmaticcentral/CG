#! /usr/bin/python
import subprocess
import time, re

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *

tcga_icgc_table_correspondence = {
"ACC_somatic_mutations" : "ACC_somatic_mutations",
"ALL_somatic_mutations" : "ALL_simple_somatic",
"BLCA_somatic_mutations": "BLCA_simple_somatic",
"BRCA_somatic_mutations": "BRCA_simple_somatic",
"CESC_somatic_mutations": "CESC_simple_somatic",
"CHOL_somatic_mutations": "CHOL_somatic_mutations",
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
"MESO_somatic_mutations": "MESO_somatic_mutations",
"OV_somatic_mutations"  : "OV_simple_somatic",
"PAAD_somatic_mutations": "PACA_simple_somatic",
"PCPG_somatic_mutations": "PCPG_somatic_mutations",
"PRAD_somatic_mutations": "PRAD_simple_somatic",
"READ_somatic_mutations": "COCA_simple_somatic",
"SARC_somatic_mutations": "SARC_simple_somatic",
"SKCM_somatic_mutations": "MELA_simple_somatic",
"STAD_somatic_mutations": "GACA_simple_somatic",
"TGCT_somatic_mutations": "TGCT_somatic_mutations",
"THCA_somatic_mutations": "THCA_simple_somatic",
"THYM_somatic_mutations": "THYM_somatic_mutations",
"UCEC_somatic_mutations": "UCEC_simple_somatic",
"UCS_somatic_mutations" : "UTCA_simple_somatic",
"UVM_somatic_mutations" : "UVM_somatic_mutations"
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
pathogenic = {'stop_lost', 'inframe_deletion', 'inframe_insertion', 'stop_gained', '5_prime_UTR_premature_start_codon_gain',
                     'start_lost', 'frameshift', 'disruptive_inframe_deletion',
                     'exon_loss', 'disruptive_inframe_insertion', 'missense',
                     'splice_acceptor', 'splice_region', 'splice_donor', 'inframe'
             }


#########################################
def quotify(something):
	if not something:
		return ""
	if type(something)==str:
		return "\'"+something+"\'"
	else:
		return str(something)


#########################################
def check_location_seen(cursor, annovar_named_field):
	location_table = "locations_chrom_%s" % annovar_named_field['chr']
	qry = "select count(*) from %s where position=%s"%(location_table, annovar_named_field['start'])
	ret = search_db(cursor,qry)
	return False if not ret or ret[0][0]==0 else True;


#########################################
def check_mutation_seen(cursor, annovar_named_field):
	location_table = "mutations_chrom_%s" % annovar_named_field['chr']
	qry = "select count(*) from %s where start_position=%s "%(location_table, annovar_named_field['start'])
	qry += "and mutated_from_allele='%s' and mutated_to_allele='%s' "%(annovar_named_field['ref'],annovar_named_field['alt'])
	ret = search_db(cursor,qry)
	return False if not ret or ret[0][0]==0 else True;



#########################################
def process_table(home, cursor, tcga_table, icgc_table, already_deposited_samples):

	# make a workdir and move there
	tumor_short = tcga_table.split("_")[0]

	qry  = "select * from tcga.%s " % tcga_table

	exit()

	# store the annotated input to icgc tables - *_simple_somatic, mutations_chrom*, and locations_chrom_*
	store_variant (cursor, icgc_table)


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
		exit()

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
