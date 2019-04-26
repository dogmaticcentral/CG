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
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#
import subprocess
import time, re

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *
from icgc_utils.icgc import *
from icgc_utils.annovar import *

from config import Config

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
pathogenic = {'stop_lost', 'inframe_deletion', 'inframe_insertion', 'stop_gained', '5_prime_UTR_premature_start_codon_gain',
			  'start_lost', 'frameshift', 'disruptive_inframe_deletion',
					 'exon_loss', 'disruptive_inframe_insertion', 'missense',
					 'splice_acceptor', 'splice_region', 'splice_donor', 'inframe', 'splice'
			 }


################################################################
def create_icgc_table(cursor, tcga_table):

	icgc_table = tcga_table.split("_")[0]+"_simple_somatic"
	if check_table_exists(cursor, 'icgc', icgc_table): return
	make_specimen_table(cursor, 'icgc', icgc_table)


#########################################
mutation_annot_pattern = re.compile('(\D+)(\-*\d+)(\D+)')
#########################################
def parse_mutation (mutation):

	match_return = re.match(mutation_annot_pattern, mutation)
	if not match_return: return [None, None, None]
	mut_from = match_return.group(1)
	mut_to   = match_return.group(3)
	mut_position = int (match_return.group(2))
	return [mut_position, mut_from, mut_to]


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
		# not sure here why some entries pop up as bytes rather than str
		# (this is python3 issue; does not appear in python)
		(tumor_sample_barcode, chromosome, start_position, end_position, reference_allele,
			tumor_seq_allele1, tumor_seq_allele2, assembly) = \
			[str(entry, 'utf-8') if type(entry)==bytes else str(entry) for entry in row]
		if tumor_sample_barcode in already_deposited_samples: continue
		differing_allele = tumor_seq_allele1
		if differing_allele==reference_allele: differing_allele = tumor_seq_allele2
		# somebody in  the TCGA decided to innovate and mark the insert with the before- and after- position
		# Annovar expects both numbers to be the same in such case
		if reference_allele=="-": end_position=start_position
		outrow = "\t".join([chromosome, start_position, end_position, reference_allele, differing_allele])
		outf[assembly].write(outrow+"\n")

	for assembly in assemblies:
		outf[assembly].close()
		subprocess.call(["bash","-c", "sort %s | uniq > %s.tmp" % (outfname[assembly], outfname[assembly])])
		os.rename(outfname[assembly]+".tmp", outfname[assembly])

	return outfname


##########
def get_icgc_mutation_type(start, end, ref, alt):
	frame_consequences = set([])
	if ref=='-':
		if len(alt)%3==0:
			frame_consequences.add('inframe')
		else:
			frame_consequences.add('frameshift')
		return 'insertion', frame_consequences
	if alt=='-':
		if len(ref)%3==0:
			frame_consequences.add('inframe')
		else:
			frame_consequences.add('frameshift')
		return 'deletion', frame_consequences
	if start==end:
		return 'single', frame_consequences
	return 'multiple', frame_consequences


#########################################
def parse_annovar_fields(cursor, avfile, annovar_named_field):

	# I should probably get rid of mutation_type columns - it is completely derivable from other columns
	# note though that it might complement consequence notation which only says inframe or frameshift
	tr_relative = []
	aa_mutation = []
	consequences = set([])
	mutation_type, frame_consequences = get_icgc_mutation_type(annovar_named_field['start'], annovar_named_field['end'],
	                                                           annovar_named_field['ref'], annovar_named_field['alt'])
	# icgc fields in locations_chrom_*: position, gene_relative, transcript relative
	# gene relative is in two annovar fields: func and gene
	gene_relative = dict(list(zip(annovar_named_field['gene'].split(";"),annovar_named_field['func'].split(";") )))
	gene_relative_string = ";".join(["%s:%s"%(k,gene_location_annovar_to_icgc(v))
	                                 for k,v in gene_relative.items()])
	# what to do with introns? can introns be  assigned to a transcript?
	# tcga thinks they can, while annovar things they are only assignable to a gene as a whole (I'd agree)
	# I am disregarding intronic mutations as zero impact
	# as an ad hoc measure, I will assign the annotation to the canonical transcript
	if 'intronic' in annovar_named_field['func']:
		consequences.add('intronic')
		for gene_stable_id, annot in gene_relative.items():
			if annot != 'intronic': continue
			canonical_transcript_id = gene_stable_id_2_canonical_transcript_id(cursor, gene_stable_id)
			tr_relative.append("{}:{}".format(canonical_transcript_id, 'intronic'))
	# transcript relative contains a bit more info than I need - perhaps if I was doing it
	# again I would stick to annovar annotation, but not now
	# I would in any case change the names in the annovar header - they do not reflect the content
	for trrel in annovar_named_field['aachange'].split(","):
		fields = trrel.split(":") # the content of each field is not fixed
		[enst,cdna_change, aa_change, annotation] = [None]*4
		for field in fields:
			if field[:4] == "ENSG":
				# just drop, we already have that info
				pass
			elif field[:4] == "ENST":
				enst = field
			elif field[:2] == "c.":
				cdna_change = field
				pass
			elif field[:2] == "p.":
				aa_change = field
			elif field[:4] == "exon":
				annotation = 'exon'
			elif field[:4] == "UTR3":
				annotation = 'UTR3'
			elif field[:4] == "UTR5":
				annotation = 'UTR5'

		if not enst:
			if aa_change:
				print("aa change without specified ENST")
				print(avfile)
				print(fields)
				exit()
			else:
				continue

		if annotation=='exon' and cdna_change and ('-' in cdna_change or '+' in cdna_change):
			annotation = 'splice_region'
			consequences.add('splice_region')

		tr_relative.append("{}:{}".format(enst, annotation))

		if aa_change:
			aa_mutation.append("{}:{}".format(enst, aa_change[2:]))
			if aa_change[-1] in ['X','x','*']:
				consequences.add('stop_gained')
			elif aa_change[-2:].lower()=='fs':
				consequences.add('frameshift')
			elif aa_change[-2:].lower() in ['ins','del']:
				pass
			elif aa_change[2]==aa_change[-1]:
				consequences.add('synonymous')
			else:
				consequences.add('missense')

	tr_relative_string  = ";".join(list(tr_relative))
	if 'exon' in tr_relative_string: consequences |= frame_consequences
	consequences_string = ";".join(list(consequences))
	aa_mutation_string = ";".join(list(aa_mutation))

	return [mutation_type, consequences_string, aa_mutation_string,
	        gene_relative_string, tr_relative_string]


#########################################
def check_location_seen(cursor, annovar_named_field, lock_alias=None):
	location_table = "locations_chrom_%s" % annovar_named_field['chr']
	qry = "select count(*) from %s " % location_table
	if lock_alias:
		qry += "as %s "% lock_alias
	qry += "where position=%s" % annovar_named_field['start']
	ret = search_db(cursor,qry)
	return False if not ret or ret[0][0]==0 else True


#########################################
def check_mutation_seen(cursor, annovar_named_field, lock_alias=None):
	mutation_table = "mutations_chrom_%s" % annovar_named_field['chr']
	qry = "select count(*) from %s " % mutation_table
	if lock_alias:
		qry += "as %s "% lock_alias #
	qry += "where start_position=%s "% annovar_named_field['start']
	qry += "and mutated_from_allele='%s' and mutated_to_allele='%s' "%(annovar_named_field['ref'],annovar_named_field['alt'])
	ret = search_db(cursor,qry)
	return False if not ret or ret[0][0]==0 else True

#########################################
def store_location(cursor, annovar_named_field, gene_relative_string, tr_relative_string):
	location_table = "locations_chrom_%s" % annovar_named_field['chr']
	# not sure if this one can go into race cond too,  but now that I am here
	# I will lock anyway
	lock_alias = "locslock"
	qry = "lock tables %s write, %s as %s read" % (location_table,location_table, lock_alias)
	search_db(cursor,qry)

	# now back to business: have we seen this mutation:
	if not check_location_seen(cursor, annovar_named_field, lock_alias):
		# MySQL  will auto-convert data types as best it can. stil let's try to convert the position here
		named_fields = {'position': int(annovar_named_field['start']),
						'gene_relative': gene_relative_string,
						'transcript_relative': tr_relative_string}

		store_without_checking(cursor, location_table, named_fields)
	# unlock
	qry = "unlock tables"
	search_db(cursor,qry)
	return


#########################################
def store_mutation(cursor, annovar_named_field, assembly, mutation_type, consequences_string, aa_change, pathogenicity_estimate):
	mutation_table = "mutations_chrom_%s" % annovar_named_field['chr']

	# I've seen this get into race condition  - I will have to lock here (? what about myISAM not deadlocking?)
	# for some reason you cannot refer to a locked table multiple times in a single query using the same name.
	#  Use aliases instead, and obtain a separate lock for the table and each alias
	# that is you do something like "select * from t as myalias"
	lock_alias = "mutslock"
	qry = "lock tables %s write, %s as %s read" % (mutation_table,mutation_table, lock_alias)
	search_db(cursor,qry)

	# now back to business: have we seen this mutation:
	if not check_mutation_seen(cursor, annovar_named_field, lock_alias):
		# find the last used id and create a new one
		qry  = "select icgc_mutation_id from %s as %s "  % (mutation_table, lock_alias)
		qry += "where icgc_mutation_id like 'MUT_%' order by icgc_mutation_id desc limit 1"
		ret = search_db(cursor,qry)
		ordinal = 1
		if ret:
			ordinal = int( ret[0][0].split("_")[-1].lstrip("0") ) + 1
		new_id = "MUT_%s_%08d"%(annovar_named_field['chr'],ordinal)
		ref_allele = annovar_named_field['ref']
		if len(ref_allele)>200: ref_allele=ref_allele[:200]+"etc"
		to_allele  = annovar_named_field['alt']
		if len(to_allele)>200: to_allele=to_allele[:200]+"etc"

		named_fields = {'icgc_mutation_id':	new_id,
						'start_position': int(annovar_named_field['start']),
						'end_position': int(annovar_named_field['end']),
						'assembly': assembly,
						'mutation_type':mutation_type,
						'reference_genome_allele': ref_allele,
						'mutated_from_allele': ref_allele,
						'mutated_to_allele': to_allele,
						'aa_mutation':aa_change,
						'consequence':consequences_string,
						'pathogenicity_estimate':pathogenicity_estimate,
						'reliability_estimate':1}
		store_without_checking(cursor, mutation_table, named_fields, verbose=False)

	# unlock
	qry = "unlock tables"
	search_db(cursor,qry)
	return


#########################################
#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./icgc_utils/kernprof.py -l 30_reorganize_tcga_mutations_and_locations.py
# followed by
# python3 -m line_profiler 30_reorganize_tcga_mutations_and_locations.py.lprof
# see here https://github.com/rkern/line_profiler#line-profiler
# the reason I am using local kernprof.py is that I don't know where pip
# installed its version (if anywhere)
#@profile
#########################################
def store_annotation (cursor, tcga_table, avoutput):
	# store location and mutation info
	# the do another sweep through the tcga table to associate with variants
	# in principle, we should be working with icgc here
	switch_to_db(cursor, 'icgc')
	for assembly, avfile in avoutput.items():
		no_lines = int(subprocess.check_output(["bash","-c", "wc -l %s"%avfile]).split()[0])
		inf = open (avfile, "r")
		header_fields = None
		ct = 0
		time0 = time.time()
		for line in inf:
			ct += 1
			if (ct%1000==0):
				print("%30s   %6d lines out of %6d  (%d%%)  %d min" % \
					  (tcga_table, ct, no_lines, float(ct)/no_lines*100, float(time.time()-time0)/60))
			if not header_fields:
				header_fields = parse_annovar_header_fields(line)
				if not header_fields:
					print("Error parsing annovar: no header line found")
					exit()
				continue
			fields = line.rstrip().split('\t')
			annovar_named_field = dict(list(zip(header_fields,fields)))
			# do I have this location already?
			location_seen = check_location_seen(cursor, annovar_named_field)
			# do I have this mutation already?
			mutation_seen = check_mutation_seen(cursor, annovar_named_field)
			# if yes to both, move on
			if location_seen and mutation_seen: continue
			ret  = parse_annovar_fields(cursor, avfile, annovar_named_field)
			[mutation_type, consequences_string, aa_change, gene_relative_string, tr_relative_string] = ret
			pathogenicity_estimate = 0
			for description in pathogenic:
				if description in consequences_string +";"+ tr_relative_string:
					pathogenicity_estimate=1
					break
			if not location_seen:
				#print("storing location:", annovar_named_field)
				store_location(cursor, annovar_named_field, gene_relative_string, tr_relative_string)
			if not mutation_seen:
				store_mutation(cursor, annovar_named_field, assembly, mutation_type,
							   consequences_string, aa_change, pathogenicity_estimate)

		inf.close()

#########################################
def process_table(home, cursor, tcga_table, already_deposited_samples):

	# make a workdir and move there
	tumor_short = tcga_table.split("_")[0]
	workdir  = tumor_short + "_annovar"
	workpath = "{}/annovar/{}".format(home,workdir)
	if not os.path.exists(workpath): os.makedirs(workpath)
	os.chdir(workpath)

	# use tcga entries to create input for annovar
	avinput  = output_annovar_input_file (cursor, 'tcga', tcga_table, already_deposited_samples)
	avoutput = {}
	for assembly, input_file in avinput.items():
		# fork annovar process
		avoutput[assembly] = run_annovar(input_file, assembly, tcga_table)
	# store the annotated input to icgc tables - *_simple_somatic, mutations_chrom*, and locations_chrom_*
	store_annotation (cursor, tcga_table, avoutput)


#########################################
def add_tcga_diff(tcga_tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	home = os.getcwd()
	for tcga_table in tcga_tables:

		icgc_table =  tcga_icgc_table_correspondence[tcga_table]

		# where in the icgc classification does this symbol belong?
		time0 = time.time()
		print()
		print("====================")
		print("processing tcga table %s, process id %d " % (tcga_table, os.getpid()))

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
		print("not deposited:", len(samples_not_deposited))
		process_table(home, cursor, tcga_table, already_deposited_samples)
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

	cursor.close()
	db.close()

	number_of_chunks = 20  # myISAM does not deadlock
	parallelize(number_of_chunks, add_tcga_diff, tcga_tables, [])

#########################################
if __name__ == '__main__':
	main()
