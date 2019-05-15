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

from icgc_utils.processes import *
from icgc_utils.tcga import *
from icgc_utils.annovar import *
from icgc_utils.common_queries import *
from config import Config

# this is set literal
pathogenic = {'stop_lost', 'inframe_deletion', 'inframe_insertion', 'stop_gained', '5_prime_UTR_premature_start_codon_gain',
			  'start_lost', 'frameshift', 'disruptive_inframe_deletion',
					 'exon_loss', 'disruptive_inframe_insertion', 'missense',
					 'splice_acceptor', 'splice_region', 'splice_donor', 'inframe', 'splice'
			 }



#########################################
def store_specimen_info(cursor, tumor_short, donor_id, tcga_barcode):
	fixed_fields  = {'icgc_donor_id':donor_id, 'icgc_specimen_id': tcga_barcode}
	update_fields = {'specimen_type':specimen_type_from_TCGA_barcode(tcga_barcode)}
	store_or_update(cursor, "icgc.%s_specimen"%tumor_short, fixed_fields, update_fields)
	return


#########################################
def store_donor_info(cursor, tumor_short, donor_id, tcga_barcode):
	fixed_fields  = {'icgc_donor_id':donor_id, 'submitted_donor_id': tcga_sample2tcga_donor(tcga_barcode)}
	update_fields = None
	store_or_update(cursor, "icgc.%s_donor"%tumor_short, fixed_fields, update_fields)
	return



#########################################
def variant_key(position_translation,  chromosome, start_position, end_position,
				reference_allele, tumor_seq_allele1, tumor_seq_allele2):
	differing_allele = tumor_seq_allele1
	if differing_allele==reference_allele: differing_allele = tumor_seq_allele2
	# somebody in  the TCGA decided to innovate and mark the insert with the before- and after- position
	# Annovar expects both numbers to be the same in such case
	if reference_allele=="-": end_position=start_position
	start_position_translated = position_translation[start_position]
	end_position_translated   = position_translation[end_position]

	return [chromosome, str(start_position_translated), str(end_position_translated), reference_allele, differing_allele]


#########################################
def output_annovar_input_file (cursor, tcga_table, chromosome, position_translation, target_samples, label):

	outfname = "%s_%s_%s.avin" % (tcga_table, chromosome, label)
	# danger zone: this file better does not exist if it is outdated
	if os.path.exists(outfname) and os.path.getsize(outfname)!=0:
		print("\t\t %s found"%outfname)
		return outfname

	meta_table_name = tcga_table.split("_")[0] + "_mutations_meta"

	# get the info that annovar needs
	qry  = "select s.tumor_sample_barcode, s.chromosome, s.start_position, s.end_position, "
	qry += "s.reference_allele, s.tumor_seq_allele1, s.tumor_seq_allele2, m.assembly  "
	qry += "from tcga.%s s, tcga.%s m " % (tcga_table, meta_table_name)
	qry += "where s.chromosome='%s' and  s.meta_info_id=m.id" % chromosome
	rows = hard_landing_search(cursor, qry)

	#target_donors = [tcga_sample2tcga_donor(s) for s in target_samples_or_donors]
	outf = open(outfname, 'w')

	for row in rows:
		# not sure here why some entries pop up as bytes rather than str
		# (this is python3 issue; does not appear in python)
		(tumor_sample_barcode, chromosome, start_position, end_position, ref_allele,
			tumor_allele1, tumor_allele2, assembly) = \
			[str(entry, 'utf-8') if type(entry)==bytes else str(entry) for entry in row]
		if not tumor_sample_barcode in target_samples: continue
		key_data = variant_key(position_translation[assembly][chromosome], chromosome, start_position, end_position,
								ref_allele, tumor_allele1, tumor_allele2)
		outrow = "\t".join(key_data)
		outf.write(outrow+"\n")

	outf.close()
	subprocess.call(["bash","-c", "sort %s | uniq > %s.tmp" % (outfname, outfname)])
	os.rename(outfname+".tmp", outfname)

	return outfname


#########################################
def intragenic(gene_relative_string):
	intra_gene_relative = set([])
	for ensid in gene_relative_string.split(";"):
		if ":" in ensid:
			ensid,loc = ensid.split(":")
			if loc=='intragenic':intra_gene_relative.add(ensid)
		else:
			intra_gene_relative.add(ensid)
	if len(intra_gene_relative)==0: return None

	return ";".join(list(intra_gene_relative))

#############
def transcript_related(tr_relative_string):
	return 'ENST' in tr_relative_string

#############
def store_location(cursor, location_table, location):

	switch_to_db(cursor, 'icgc')
	[start, end,  gene_relative_string, tr_relative_string] = location
	# decision Apr 2019: store only intragenic
	intra_gene_relative = intragenic(gene_relative_string)
	transcript_relative = transcript_related(tr_relative_string)
	if not intra_gene_relative: return "not intragenic"
	if not transcript_relative: return "not transcript-related"
	named_fields = {'position': start,
					'gene_relative': gene_relative_string,
					'transcript_relative': tr_relative_string}
	# position is the key in location table - ignore if exists
	store_without_checking(cursor, location_table, named_fields, database='icgc', ignore=True)
	# TODO: do this some other way - the start and end position might not be in the same gene-relative location
	if end !=start:
		named_fields['position']=end
		store_without_checking(cursor, location_table, named_fields, database='icgc', ignore=True)
	return "ok"

#########################################
def create_new_mutation_id(cursor, chromosome, mutation_table):
	# find the last used id and create a new one
	qry  = "select icgc_mutation_id from %s "  % mutation_table
	qry += "where icgc_mutation_id like 'MUT_%' order by icgc_mutation_id desc limit 1"
	ret = search_db(cursor,qry)
	ordinal = 1
	if ret:
		ordinal = int( ret[0][0].split("_")[-1].lstrip("0") ) + 1
	return  "MUT_%s_%08d"%(chromosome,ordinal)

def make_pathogenicity_estimate(consequences_string, tr_relative_string):
	pathogenicity_estimate=0
	for description in pathogenic:
		if description in consequences_string +";"+ tr_relative_string:
			pathogenicity_estimate=1
			break
	return pathogenicity_estimate

###########
def store_mutation(cursor, mutation_table, var_key, mutation):
	[chromosome,  start, end , ref, alt] = var_key.split("_")
	[mutation_type, consequences_string, aa_mutation_string, gene_relative_string, tr_relative_string] = mutation

	# check mutation seen
	qry  = "select icgc_mutation_id, pathogenicity_estimate from %s " % mutation_table
	qry += "where start_position=%s " % start
	qry += "and mutated_from_allele='%s' and mutated_to_allele='%s' "%(ref, alt)
	qry += "limit 1" # if we have stored this mutation more than once, this is not a place to deal with it
	ret = error_intolerant_search(cursor,qry)
	if ret and ret[0]:
		#print("mutation found, ", ret[0])
		return ret[0]

	# if we have not seen this mutation yet, create new mutation id and store
	new_id = create_new_mutation_id(cursor, chromosome, mutation_table)
	pathogenicity_estimate = make_pathogenicity_estimate(consequences_string, tr_relative_string)
	named_fields = {'icgc_mutation_id':	new_id,
					'start_position': int(start),
					'end_position': int(end),
					'assembly': Config.ref_assembly,
					'mutation_type':mutation_type,
					'reference_genome_allele': ref,
					'mutated_from_allele': ref,
					'mutated_to_allele': alt,
					'aa_mutation':aa_mutation_string,
					'consequence':consequences_string,
					'pathogenicity_estimate':pathogenicity_estimate,
					'reliability_estimate':1}
	#print("storing without checking", mutation_table, named_fields)
	store_without_checking(cursor, mutation_table, named_fields, verbose=False)
	return [new_id, pathogenicity_estimate]

#########################################
def store_variant(cursor, icgc_variant_table, donor_id, specimen_id, sample_barcode, var_key, mutation_id, pathogenicity_estimate):
	[chromosome,  start, end , ref, alt] = var_key.split("_")
	#have we stored this by any chance?
	qry  = "select count(*) from icgc.%s " % icgc_variant_table
	qry += "where icgc_mutation_id='%s' " % mutation_id
	qry += "and submitted_sample_id='%s'" %  sample_barcode
	if error_intolerant_search(cursor,qry)[0][0]>0:
		# print("variant stored already")
		# print(qry)
		# print(sample_barcode, var_key, mutation_id, pathogenicity_estimate)
		#exit()
		return
	# if variants from this same donor exist under different sample/specimen heading we will have to
	# reove duplicates downstream

	# fill store hash
	store_fields = {
		'icgc_mutation_id': mutation_id,
		'chromosome': chromosome,
		'icgc_donor_id': donor_id,
		'icgc_specimen_id': specimen_id,
		'icgc_sample_id': sample_barcode,
		'submitted_sample_id': sample_barcode,
		'tumor_genotype': "{}/{}".format(ref,alt),
		'pathogenicity_estimate': pathogenicity_estimate,
		'reliability_estimate': 1
	}
	# store - we actually checked above
	#print("storing without checking", store_fields)
	store_without_checking(cursor, icgc_variant_table, store_fields, verbose=False, database='icgc')

	return

#########################################
def get_donor_id(cursor, icgc_table, sample_id):
	qry  = "select distinct(icgc_donor_id) from icgc.%s " % icgc_table
	qry += "where submitted_sample_id='%s'" % sample_id
	ret = error_intolerant_search(cursor, qry)
	if ret:
		donor_ids = [r[0] for r in ret]
		#print(donor_ids) # wht do I do in this case?
		donor_id = donor_ids[0]
	else:
		# use tcga id
		donor_id = tcga_sample2tcga_donor(sample_id)
		#print("using tcga id", donor_id)
	return donor_id

#########################################
def get_specimen_id(cursor, icgc_table, sample_id):
	qry  = "select distinct(icgc_specimen_id) from icgc.%s " % icgc_table
	qry += "where submitted_sample_id='%s'" % sample_id
	ret = error_intolerant_search(cursor, qry)
	if ret:
		specimen_ids = [r[0] for r in ret]
		#print(specimen_ids) # wht do I do in this case?
		specimen_id = specimen_ids[0]
	else:
		# use tcga id
		specimen_id = tcga_sample2tcga_donor(sample_id)
		#print("using tcga id", specimen_id)
	return specimen_id

#########################################
def resolve_ids(cursor, icgc_table, tumor_sample_barcode, flag):
	donor_id = None
	specimen_id  = None
	if flag=='donors':
		donor_id = tcga_sample2tcga_donor(tumor_sample_barcode)
		specimen_id = specimen_id_from_TCGA_barcode(tumor_sample_barcode)
	elif flag=='vars':
		donor_id = get_donor_id(cursor, icgc_table, tumor_sample_barcode)
		specimen_id = get_specimen_id(cursor, icgc_table, tumor_sample_barcode)
	else:
		print("unrecognized flag", flag)
		exit()
	return donor_id, specimen_id

#########################################
def process_samples(cursor, tcga_table, icgc_table, chromosome, target_samples, position_translation, flag):

	print("\t", flag, tcga_table, chromosome)
	avinput  = output_annovar_input_file(cursor, tcga_table, chromosome,
										position_translation, target_samples, flag)
	avoutput = run_annovar(avinput, Config.ref_assembly, tcga_table)
	annovar_annotation = parse_annovar(cursor, avoutput)
	tumor_short = icgc_table.split("_")[0]
	meta_table_name = tcga_table.split("_")[0] + "_mutations_meta"
	for tumor_sample_barcode in target_samples:
		#print(tumor_sample_barcode)
		donor_id, specimen_id = resolve_ids(cursor, icgc_table, tumor_sample_barcode, flag)
		qry  = "select  s.chromosome, s.start_position, s.end_position, "
		qry += "s.reference_allele, s.tumor_seq_allele1, s.tumor_seq_allele2, m.assembly  "
		qry += "from tcga.%s s, tcga.%s m " % (tcga_table, meta_table_name)
		qry += "where s.chromosome='%s' " % chromosome
		qry += "and  s.tumor_sample_barcode = '%s' " %  tumor_sample_barcode
		qry += "and  s.meta_info_id=m.id"
		ret = error_intolerant_search(cursor, qry)

		if not ret: continue
		for columns in  ret:
			clean_columns = [str(entry, 'utf-8') if type(entry)==bytes else str(entry) for entry in columns]
			(chromosome, start_position, end_position, ref_allele,
					tumor_allele1, tumor_allele2, assembly) = clean_columns
			if flag=='donors':
				# this is redundant, but then the specimen table is small ...
				store_specimen_info(cursor, tumor_short, donor_id, tumor_sample_barcode)
				store_donor_info(cursor, tumor_short, donor_id, tumor_sample_barcode)

			# (remind me to get rid of annovar)
			# unpacking  the list into positional arguments using asterisk
			key_data = variant_key(position_translation[assembly][chromosome], *(clean_columns[:-1]))
			dict_key = "_".join(key_data)
			if not dict_key in annovar_annotation:
				print(dict_key, "not found in annovar_annotation. Is the annovar dir up to date?")
				exit()
			ret_flag = store_location(cursor, "locations_chrom_%s"%chromosome, annovar_annotation[dict_key]['loc'])
			if ret_flag != "ok":
				#print(ret_flag)
				continue
			[mutation_id, pathogenicity_estimate] = \
				store_mutation(cursor,"mutations_chrom_%s"%chromosome, dict_key, annovar_annotation[dict_key]['mut'])
			store_variant(cursor, icgc_table, donor_id, specimen_id, tumor_sample_barcode, dict_key, mutation_id, pathogenicity_estimate)
	return


#########################################
def locate_complementary_ids(cursor, tcga_table, icgc_table):

	print(tcga_table, icgc_table)

	# sample ids in TCGA
	qry  = "select distinct(tumor_sample_barcode) from tcga.%s " % tcga_table
	ret =  error_intolerant_search(cursor, qry)
	if not ret: return
	tcga_sample_ids_in_tcga = set([r[0] for r in ret])
	tcga_donor_ids_in_tcga = set([tcga_sample2tcga_donor(s) for s in tcga_sample_ids_in_tcga])

	# find submitted sample ids in ICGC
	qry  = "select distinct(submitted_sample_id) from icgc.%s " % icgc_table
	qry += "where (submitted_sample_id like 'TCGA%' or submitted_sample_id like 'TARGET%')"
	ret  =  error_intolerant_search(cursor, qry)
	if not ret:
		tcga_sample_ids_in_icgc = ()
		tcga_donor_ids_in_icgc = ()
	else:
		tcga_sample_ids_in_icgc = set([r[0] for r in ret])
		tcga_donor_ids_in_icgc  = set([tcga_sample2tcga_donor(s) for s in tcga_sample_ids_in_icgc])

	donor_ids_not_in_icgc = tcga_donor_ids_in_tcga.difference(tcga_donor_ids_in_icgc)
	sample_ids_already_in_icgc = tcga_sample_ids_in_tcga.intersection(tcga_sample_ids_in_icgc)

	# I need this to be ordered, at least while debugging:
	donor_ids_not_in_icgc = sorted(list(donor_ids_not_in_icgc))
	sample_ids_already_in_icgc = sorted(list(sample_ids_already_in_icgc))

	return [donor_ids_not_in_icgc, sample_ids_already_in_icgc]

####################
def find_preferred_samples_for_donor(cursor, tcga_table, tcga_donor_id):
	qry  = "select distinct(tumor_sample_barcode) from tcga.%s " % tcga_table
	qry += "where tumor_sample_barcode like '%s%%'" % tcga_donor_id
	barcodes = [r[0] for r in hard_landing_search(cursor,qry)]
	barcodes_by_type = {}
	for barcode in barcodes:
		spec_type =  specimen_type_from_TCGA_barcode(barcode)
		if not spec_type in barcodes_by_type: barcodes_by_type[spec_type] = []
		barcodes_by_type[spec_type].append(barcode)

	# possible types 'Primary' 'Recurrent' 'Metastatic' 'Normal'  'Cell_line'
	# the order of preference; we do not want normal in any case
	for spec_type in ['Primary', 'Recurrent', 'Metastatic', 'Cell_line']:
		if spec_type in barcodes_by_type:
			return barcodes_by_type[spec_type]

	return None

#####################
def find_preferred_samples(cursor, tcga_table, donor_ids_not_in_icgc):

	preferred_samples = []
	for tcga_donor_id in donor_ids_not_in_icgc:
		# ideally, we would like the primary tumor sample; in any case we do not want the normal tissue
		tumor_sample_barcodes = find_preferred_samples_for_donor(cursor, tcga_table, tcga_donor_id)
		if not tumor_sample_barcodes: continue
		preferred_samples.extend(tumor_sample_barcodes)
	return preferred_samples

#########################################
def process_tables(tcga_tables, other_args):
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	home   = os.getcwd()
	for tcga_table in tcga_tables:
		# make a workdir and move there
		tumor_short = tcga_table.split("_")[0]
		workdir  = tumor_short + "_annovar"
		workpath = "{}/annovar/{}".format(home,workdir)
		if not os.path.exists(workpath): os.makedirs(workpath)
		os.chdir(workpath)
		# translate all positions en masse
		position_translation = get_position_translation(cursor, tcga_table, Config.ref_assembly)

		# what is the icgc table we eill be storing to?
		icgc_table = Config.tcga_icgc_table_correspondence[tcga_table]
		if not icgc_table: icgc_table = tcga_table.replace("somatic_mutations", "simple_somatic")

		# TCGA has something called sample barcode, which containus info about the donor,
		# as well as the sample (and implicitly, the specimen)
		#     There are two cases we are interested in - either we have not seen this donor at all,
		# or we we have the sample, but not all variants. The reason fo the latter is that
		# in this pipeline TCGA came first, and we needed to have all the variants we have considered in the round one.
		# Were they dropped for a reason? The ICGC is not documented well enough to know, and TCGA
		# largely does not have the sequencing depth info.
		#     If we have the donor, but different sample, we will pass, assuming that ICGC know what they are doing,
		# and the sample they have analyzed is somehow superior to what is deposited in TCGA.
		[donor_ids_not_in_icgc, sample_ids_already_in_icgc] = locate_complementary_ids(cursor, tcga_table, icgc_table)
		preferred_samples_for_compl_donors = find_preferred_samples(cursor, tcga_table, donor_ids_not_in_icgc)
		# process info on -per chromosome basis
		chromosomes = [str(i) for i in range(1,23)] + ['X','Y']
		for chrom in chromosomes:
			process_samples(cursor, tcga_table, icgc_table, chrom, preferred_samples_for_compl_donors, position_translation, 'donors')
			#process_samples(cursor, tcga_table, icgc_table, chrom, sample_ids_already_in_icgc, position_translation, 'vars')


	cursor.close()
	db.close()

	return


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

	number_of_chunks = 1
	tcga_tables = ['THCA_somatic_mutations']
	parallelize(number_of_chunks, process_tables, tcga_tables, [])


#########################################
if __name__ == '__main__':
	main()

