#! /usr/bin/python3
#
# This source code is part of icgc, an ICGC processing pipeline.
# 
# Icgc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Icgc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#

import time

from icgc_utils.common_queries import *
from icgc_utils.processes import *
from icgc_utils.icgc import *
from config import Config

tcga_icgc_table_correspondence = {
	"ACC_somatic_mutations" :  "ACC_simple_somatic",
	"ALL_somatic_mutations" :  "ALL_simple_somatic",
	"BLCA_somatic_mutations": "BLCA_simple_somatic",
	"BRCA_somatic_mutations": "BRCA_simple_somatic",
	"CESC_somatic_mutations": "CESC_simple_somatic",
	"CHOL_somatic_mutations": "CHOL_simple_somatic",
	"COAD_somatic_mutations": "COCA_simple_somatic",
	"DLBC_somatic_mutations": "DLBC_simple_somatic",
	"ESCA_somatic_mutations": "ESAD_simple_somatic",
	"GBM_somatic_mutations" :  "GBM_simple_somatic",
	"HNSC_somatic_mutations": "HNSC_simple_somatic",
	"KICH_somatic_mutations": "KICH_simple_somatic",
	"KIRC_somatic_mutations": "KIRC_simple_somatic",
	"KIRP_somatic_mutations": "KIRP_simple_somatic",
	"LAML_somatic_mutations":  "AML_simple_somatic",
	"LGG_somatic_mutations" :  "LGG_simple_somatic",
	"LIHC_somatic_mutations": "LICA_simple_somatic",
	"LUAD_somatic_mutations": "LUAD_simple_somatic",
	"LUSC_somatic_mutations": "LUSC_simple_somatic",
	"MESO_somatic_mutations": "MESO_simple_somatic",
	"OV_somatic_mutations"  :   "OV_simple_somatic",
	"PAAD_somatic_mutations": "PACA_simple_somatic",
	"PCPG_somatic_mutations": "PCPG_simple_somatic",
	"PRAD_somatic_mutations": "PRAD_simple_somatic",
	"READ_somatic_mutations": "COCA_simple_somatic",
	"SARC_somatic_mutations": "SARC_simple_somatic",
	"SKCM_somatic_mutations": "MELA_simple_somatic",
	"STAD_somatic_mutations": "GACA_simple_somatic",
	"TGCT_somatic_mutations": "TGCT_simple_somatic",
	"THCA_somatic_mutations": "THCA_simple_somatic",
	"THYM_somatic_mutations": "THYM_simple_somatic",
	"UCEC_somatic_mutations": "UCEC_simple_somatic",
	"UCS_somatic_mutations" : "UTCA_simple_somatic",
	"UVM_somatic_mutations" :  "UVM_simple_somatic"
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
def check_location_stored(cursor, tcga_named_field):
	location_table = "icgc.locations_chrom_%s" % tcga_named_field['chromosome']
	qry = "select count(*) from %s where position=%s"%(location_table, tcga_named_field['start_position'])
	ret = search_db(cursor,qry)
	if ret and len(ret)>1:
		print("problem: non-unique location id")
		print(qry)
		print(ret)
		exit()
	return False if not ret else ret[0][0]


#########################################
def find_mutation_id(cursor, tcga_named_field):
	mutation_table = "mutations_chrom_%s" % tcga_named_field['chromosome']
	qry = "select icgc_mutation_id, pathogenic_estimate from icgc.%s where start_position=%s "%(mutation_table, tcga_named_field['start_position'])
	reference_allele = tcga_named_field['reference_allele']
	differing_allele = tcga_named_field['tumor_seq_allele1']
	if differing_allele == reference_allele: differing_allele = tcga_named_field['tumor_seq_allele2']
	if len(reference_allele)>200: reference_allele=reference_allele[:200]+"etc"
	if len(differing_allele)>200: differing_allele=differing_allele[:200]+"etc"
	qry += "and mutated_from_allele='%s' and mutated_to_allele='%s' "%(reference_allele,differing_allele)
	ret = search_db(cursor,qry)

	if not ret:
		print("problem: no return for")
		print(qry)
		exit() # brilliant idea in case of multithreading
	if len(ret)>1:
		print("problem: non-unique mutation id")
		print(qry)
		print(ret)
		exit()
	return False if not ret else ret[0]


#########################################
def construct_id(cursor, icgc_table, lock_alias, lock_alias_2,  tcga_donor_id):

	tumor_short = icgc_table.split("_")[0]
	# construct donor id
	qry  = "select icgc_donor_id from icgc.%s as %s "  % (icgc_table, lock_alias)
	qry += "where icgc_donor_id like 'DOT_%' order by icgc_donor_id desc limit 1"
	ret = search_db(cursor,qry)
	if ret and ret[0] and type(ret[0][0])==str and 'error' in ret[0][0].lower():
		search_db(cursor,qry, verbose=True)
		exit(1)
	ordinal = 1
	if ret:
		ordinal = int( ret[0][0].split("_")[-1].lstrip("0") ) + 1
	new_donor_id = "DOT_%s_%05d"%(tumor_short,ordinal)

	last_id = 0
	# locked read takes alias, locked write (store_fields below) does not
	qry = "select id from icgc.{}_donor as {} order by id desc limit 1".format(tumor_short, lock_alias_2)
	ret = search_db(cursor,qry)
	if ret:
		if type(ret[0][0])==str and 'Error' in ret[0][0]:
			search_db(cursor,qry, verbose="True")
			exit()
		last_id = ret[0][0]

	store_fields = {"id":(last_id+1), "icgc_donor_id":new_donor_id, "submitted_donor_id": tcga_donor_id}
	store_without_checking(cursor, "%s_donor"%tumor_short, store_fields, verbose=False, database='icgc')

	return new_donor_id


#########################################
def tcga_sample2donor(s):
	return "-".join(s.split("-")[:3])


#########################################
def store_variant(cursor,  tcga_named_field, mutation_id, pathogenic_estimate, icgc_table, id_resolution):


	# thread parallelization goes over tcga tables - no guarantee there won't be race condition
	# for icgc tables if two tcga's map to the same icgc
	# lock table
	tumor_short = icgc_table.split("_")[0]
	lock_alias = "{}varslock".format(tumor_short)
	lock_alias_2 = "{}donorslock".format(tumor_short)
	# You cannot refer to a locked table multiple times in a single query using the same name. Use aliases instead.
	# But then not that once  you lock a table using an alias, you must refer to it in your statements using that alias.
	# So her for example wrtie take no alias, but read does
	qry = "lock tables icgc.%s write, icgc.%s as %s read. " % (icgc_table, icgc_table, lock_alias)
	qry +="icgc.%s_donor write,  icgc.%s_donor as %s read" % (tumor_short, tumor_short, lock_alias_2)
	search_db(cursor,qry)

	#have we stored this by any chance?
	qry  = "select submitted_donor_id from icgc.%s " % icgc_table
	qry += "where icgc_mutation_id='%s' " % mutation_id
	ret = search_db(cursor,qry)
	if ret and (tcga_named_field['tumor_sample_barcode'] in [r[0] for r in ret]):
		#print "variant found"
		pass
	else:
		tcga_participant_id = tcga_sample2donor(tcga_named_field['tumor_sample_barcode'])
		if tcga_participant_id in id_resolution:
			new_donor_id = id_resolution[tcga_participant_id]
		else:
			new_donor_id = construct_id(cursor, icgc_table, lock_alias, lock_alias_2, tcga_participant_id)
			id_resolution[tcga_participant_id] = new_donor_id

		# tcga could not agree with itself in which column to place the cancer allele
		reference_allele = tcga_named_field['reference_allele']
		differing_allele = tcga_named_field['tumor_seq_allele1']
		if differing_allele == reference_allele: differing_allele = tcga_named_field['tumor_seq_allele2']
		if len(reference_allele)>200: reference_allele=reference_allele[:200]+"etc"
		if len(differing_allele)>200: differing_allele=differing_allele[:200]+"etc"

		# fill store hash
		store_fields = {
			'icgc_mutation_id': mutation_id,
			'chromosome': tcga_named_field['chromosome'],
			'icgc_donor_id': new_donor_id,
			'submitted_sample_id':tcga_named_field['tumor_sample_barcode'],
			'tumor_genotype': "{}/{}".format(reference_allele,differing_allele),
			'pathogenic_estimate': pathogenic_estimate,
			'reliability_estimate': 1,
		}

		# store
		store_without_checking(cursor, icgc_table, store_fields, verbose=False, database='icgc')

	# unlock
	qry = "unlock tables"
	search_db(cursor,qry)

	return


#########################################
def process_tcga_table(cursor, tcga_table, icgc_table, donors_already_deposited):

	standard_chromosomes = [str(i) for i in range(23)] +['X','Y']

	# make a workdir and move there
	no_rows = search_db(cursor,"select count(*) from tcga.%s"% tcga_table)[0][0]

	column_names = get_column_names(cursor,'tcga',tcga_table)
	qry = "select * from tcga.%s " % tcga_table
	ct = 0
	time0 = time.time()
	id_resolution = {}

	for row in search_db(cursor,qry):
		ct += 1
		if (ct%10000==0):
			print("%30s   %6d lines out of %6d  (%d%%)  %d min" % \
				(tcga_table, ct, no_rows, float(ct)/no_rows*100, float(time.time()-time0)/60))
		named_field = dict(list(zip(column_names,row)))

		if tcga_sample2donor(named_field['tumor_sample_barcode']) in donors_already_deposited: continue
		if not named_field['chromosome'] in standard_chromosomes: continue
		mutation_id, pathogenic_estimate = find_mutation_id(cursor, named_field)
		location_stored = check_location_stored(cursor, named_field)
		if not mutation_id or not location_stored:
			print("mutation id:", mutation_id, "location stored:", location_stored)
			#exit()
		# all clear - store
		store_variant(cursor, named_field, mutation_id, pathogenic_estimate, icgc_table, id_resolution)


#########################################
def add_tcga_diff(tcga_tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	for tcga_table in tcga_tables:

		icgc_table =  tcga_icgc_table_correspondence[tcga_table]

		# where in the icgc classification does this symbol belong?
		time0 = time.time()
		print("\n"+"-"*50+"\npid {} processing tcga table {} - will be stored in {}".\
				format(os.getpid(), tcga_table,  icgc_table))

		#tcga donors already deposited in icgc
		tumor = icgc_table.split("_")[0]
		deposited_tcga_donor_ids = []
		if check_table_exists(cursor, "icgc", "{}_donor".format(tumor)):
			qry = "select submitted_donor_id from icgc.{}_donor ".format(tumor)
			qry += "where submitted_donor_id like '{}' or  submitted_donor_id like '{}' ".format("TCGA_%", "TARGET_%")
			ret = search_db(cursor,qry, verbose=True)
			if ret: deposited_tcga_donor_ids = [r[0] for r in  ret]
			print("\t already deposited: ", len(deposited_tcga_donor_ids))
		process_tcga_table(cursor, tcga_table, icgc_table, deposited_tcga_donor_ids)
		print("\t overall time for %s: %.3f mins; pid: %d" % (tcga_table, float(time.time()-time0)/60, os.getpid()))

	cursor.close()
	db.close()

	return


#########################################
#########################################
def main():

	# divide by cancer types, because I have duplicates within each cancer type
	# that I'll resolve as I go, but I do not want the threads competing)
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='tcga' and table_name like '%_somatic_mutations'"
	tcga_tables = [field[0] for field in search_db(cursor,qry)]


	number_of_chunks = 1 # myISAM does not deadlock
	parallelize(number_of_chunks, add_tcga_diff, tcga_tables, [])


#########################################
if __name__ == '__main__':
	main()
