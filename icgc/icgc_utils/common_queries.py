
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

from icgc_utils.tcga import *
#
#########################################

def check_gene_hgnc(cursor, gene):
	qry = "select * from icgc.hgnc where approved_symbol='%s' limit 1" %gene
	ret = error_intolerant_search(cursor, qry)
	return ret!=None

def get_somatic_variant_tables(cursor):
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in error_intolerant_search(cursor,qry)]
	return tables

def gnomad_mutations (cursor, gene_symbol):

	mutations = set([])

	chromosome = find_chromosome(cursor, gene_symbol)
	#column_names
	colnames = get_column_names(cursor,"gnomad","gnomad_freqs_chr_1")

	# brute force approach seems to be fast enough for a single gene
	qry = "select * from gnomad.gnomad_freqs_chr_{} where consequences like '%|{}|%' ".format(chromosome, gene_symbol)
	ret = error_intolerant_search(cursor,qry)
	if not ret:
		print("nothing found for {}, chromosome {}".format(gene_symbol, chromosome))
		exit()
	for line in ret:
		named_fields = dict(list(zip(colnames,line)))
		if int(named_fields['variant_count'])<2: continue
		relevant_consequences = []
		for description in named_fields['consequences'].split(","):
			if not  gene_symbol in description: continue
			if not 'missense' in description: continue
			description_field = description.split("|")
			# I don't have ensembl info here - in a more through implementation one should
			# at least go for annotator here
			# for now just hope that the uniprot is canonical
			if len(description_field[2])==0: continue
			relevant_consequences.append(description)
		if len(relevant_consequences)==0: continue
		for description in relevant_consequences:
			description_field = description.split("|")
			mutation = description_field[8].split("/")[0] + description_field[7] + description_field[8].split("/")[1]
			mutations.add(mutation)

	return list(mutations)


#########################################
def count_entries(cursor, somatic_table, icgc_specimen_id):
	qry = "select count(*) from {} where icgc_specimen_id='{}' ".format(somatic_table, icgc_specimen_id)
	ret = search_db(cursor, qry)
	if not ret or ret[0][0] == 0: return 0
	return ret[0][0]


#########################################
def find_spec_id_with_max_entries(spec_ids_w_description, entries_per_specimen):
	max_count = 0
	max_spec_id = None
	for psi in spec_ids_w_description:
		[spec_id, description] = psi.split(":")
		print("\t\t ", spec_id, entries_per_specimen[spec_id], description)
		if max_count < entries_per_specimen[spec_id]:
			max_count = entries_per_specimen[spec_id]
			max_spec_id = spec_id
	return max_spec_id

#########################################
def resolve_duplicate_specimens(cursor, somatic_table, specimen_ids):
	tumor = somatic_table.split("_")[0]

	for icgc_donor_id, donor_spec_ids in specimen_ids.items():
		print("specimen ids for %s" % icgc_donor_id, donor_spec_ids)

		qry = "select icgc_specimen_id, count(*) as c  from %s " % somatic_table
		qry += "where icgc_donor_id = '%s' and reliability_estimate=1 " % icgc_donor_id
		qry += "group by  icgc_specimen_id"
		ret2 = search_db(cursor, qry)
		if not ret2:
			search_db(cursor, qry, verbose=True)
			exit(1)
		entries_per_specimen = dict(ret2)
		print(entries_per_specimen)

		qry = "select icgc_specimen_id, icgc_donor_id, specimen_type "
		qry += "from %s_specimen where icgc_specimen_id in (%s)" % \
				(tumor, ",".join(["'%s'" % id for id in entries_per_specimen.keys()]))
		ret3 = search_db(cursor, qry)
		if not ret3:
			search_db(cursor, qry, verbose=True)
			exit(1)
		other_spec_ids = set()
		primary_spec_ids = set()
		normal_spec_ids = set()
		metastatic_spec_ids = set()
		removable_ids = set(donor_spec_ids)
		for line in ret3:
			[icgc_specimen_id, icgc_donor_id, specimen_type] = line
			# here I am counting on the capitalization
			# there might say things like "Normal - tissue adjacent to primary"
			# not sure how to handle this in a completely general way
			if 'Primary' in specimen_type:
				primary_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))
			elif 'Normal' in specimen_type:
				normal_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))
			elif 'Metastatic' in specimen_type:
				metastatic_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))
			else:
				other_spec_ids.add("{}: {}".format(icgc_specimen_id, specimen_type))

		#################################
		keep_id = None
		if len(primary_spec_ids) == 0:
			print("\t no primary ids", icgc_donor_id)
			print("\t other spec ids", other_spec_ids)
			if len(metastatic_spec_ids) == 1:
				metastatic_spec_id = metastatic_spec_ids.pop().split(":")[0]
				print("there is only one metastastic tumor id: ", metastatic_spec_id)
				keep_id = metastatic_spec_id

			elif len(metastatic_spec_ids) > 1:
				max_spec_id = find_spec_id_with_max_entries(metastatic_spec_ids, entries_per_specimen)
				if max_spec_id:
					keep_id = max_spec_id
					print("using  metastastic tumor id: ", max_spec_id)
			else:  # if all else fails, use whatever specimens are available
				max_spec_id = find_spec_id_with_max_entries(other_spec_ids, entries_per_specimen)
				if max_spec_id:
					keep_id = max_spec_id
					print("using  'other' tumor id: ", max_spec_id)

		elif len(primary_spec_ids) > 1:
			print("\t multiple primary ids", icgc_donor_id)
			max_spec_id = find_spec_id_with_max_entries(primary_spec_ids, entries_per_specimen)
			if max_spec_id:
				keep_id = max_spec_id
		else:
			primary_spec_id = primary_spec_ids.pop().split(":")[0]
			print("there is only one primary id: ", primary_spec_id)
			keep_id = primary_spec_id

		##
		if keep_id:
			removable_ids.remove(keep_id)
			removable_ids_string = ",".join(["'%s'" % rid for rid in removable_ids])
			qry = "delete from {} where icgc_specimen_id in ({})".format(somatic_table, removable_ids_string)
			print(qry)
			search_db(cursor, qry)


#########################################
def resolve_duplicate_mutations(cursor, table,  duplicate_lines, verbose=True):

	colnames = get_column_names(cursor, "icgc", table)

	[max_depth, max_allele_depth] = [-1,-1]
	[max_id, max_allele_id] = [-1,-1]
	all_ids = []
	genotype = []
	path_estimate = []
	for line2 in duplicate_lines:
		if verbose: print(line2)
		named_field = dict(list(zip(colnames,line2)))
		all_ids.append( named_field["id"])
		genotype.append(named_field["tumor_genotype"])
		path_estimate.append(named_field["pathogenicity_estimate"])
		# first see what is the best total read count that we have
		if named_field["total_read_count"] and max_depth<named_field["total_read_count"]:
			max_depth = named_field["total_read_count"]
			max_id    = named_field["id"]
		# as the second line tiebreaker -- this could only happen if the total read count is null
		if named_field["mutant_allele_read_count"] and max_allele_depth<named_field["mutant_allele_read_count"]:
			max_allele_depth = named_field["mutant_allele_read_count"]
			max_allele_id    = named_field["id"]

	# the first choice: the entry with the greatest sequencing depth
	if max_id>=0 or max_allele_id>=0:
		other_ids = set(all_ids)
		if max_id>=0: # ids start from 1
			other_ids.remove(max_id)
			if verbose: print("max depth %d found at %d, removing"%( max_depth, max_id), other_ids )
		elif max_allele_id>=0:
			other_ids.remove(max_allele_id)
			if verbose: print("max allele depth %d found at %d, removing"%(max_allele_depth, max_allele_id),other_ids )
		qry = "delete from %s where id in (%s)" % (table, ",".join([str(other_id) for other_id in other_ids]))

	# we do not have the info about the depth of the sequencing,
	# but genotypes are actually the same, so it does not matter
	elif len(duplicate_lines)==2 and genotype[0]==genotype[1][::-1]: # hack to reverse a string
		#print("tumor genotypes same", genotype)
		qry = "delete from %s where id = %d" % (table, all_ids[1])

	# I am losing patience a bit here, I guess
	# if none of the entries is pathogenic (and we are takina a rather generous
	# definition of pathogenicity: frameshift, any missense, splice)
	# then just delete them all
	elif set(path_estimate)=={0}:
		#print("all entries for %s estimated irrelevant (in terms of pathogenicity)" % mega_id)
		qry = "delete from %s where id in (%s) " % (table, ",".join([str(s) for s in all_ids]))

	else: # I really cannot decide; therefore merge the annotations into the first entry and delete the rest
		qry1 = "update %s set tumor_genotype='%s' where id=%d" % (table, ";".join(genotype), all_ids[0])
		search_db(cursor,qry1)
		qry = "delete from %s where id in (%s) " % (table, ",".join([str(s) for s in all_ids[1:]]))

	search_db(cursor,qry, verbose=verbose)


#########################################
def transcript_location_cleanup(cursor, loc, gene_stable_id):
	if not loc: return ""
	if loc== "": return loc
	location = {};
	for enst_loc in loc.split(";"):
		[e, c] = enst_loc.split(":")
		location[e] = c
	enst_canonical = list_of_transcript_ids_2_canonical_transcript_id(cursor, list(location.keys()))
	if not enst_canonical: return loc
	return location[enst_canonical]


#########################################
def annotation_to_dict(aa_change):
	change = {}
	for enst_change in aa_change.split(";"):
		[e, c] = enst_change.split(":")
		change[e] = c
	return change


#########################################
def consequence_cleanup(canonical_transcript, aa_change):
	if not aa_change: return ""
	if aa_change=="": return aa_change
	if not ":" in aa_change: return aa_change
	change = annotation_to_dict(aa_change)
	if not canonical_transcript in change: return aa_change
	return change[canonical_transcript]


#########################################
def find_background_status(cursor, tumor_short, specimen, bg_gene_symbol, bg_canonical_transcript):

	# which chromosome is our background gene on?
	chromosome = find_chromosome(cursor, bg_gene_symbol)

	# g = gene
	# v = variant
	# m = mutation
	# l = location
	qry  = "select g.icgc_mutation_id, v.pathogenicity_estimate, m.consequence, m.aa_mutation, l.transcript_relative "
	qry += "from mutation2gene g, %s_simple_somatic v,  " % (tumor_short)
	qry += "mutations_chrom_%s m , locations_chrom_%s l " % (chromosome, chromosome)
	qry += "where g.gene_symbol='%s' " % bg_gene_symbol
	qry += "and v.icgc_specimen_id = '%s' " % specimen
	qry += "and g.icgc_mutation_id = v.icgc_mutation_id "
	qry += "and g.icgc_mutation_id = m.icgc_mutation_id "
	qry += "and m.start_position = l.position "
	ret = search_db(cursor,qry)
	if not ret: return ["wt",""]

	impact_estimate = "benign"
	cons = []
	for line in ret:
		[icgc_mutation_id, pathogenicity, consequence, aa_mutation, transcript_relative] = line
		if pathogenicity==1:  impact_estimate = "pathogenic"
		if transcript_relative!=None and 'splice' in transcript_relative:
			clean = consequence_cleanup(cursor,transcript_relative)
			if clean and ":" in clean: clean = clean.split(":")[1]
			cons.append(clean)
			continue
		if consequence and aa_mutation:
			aa_change = consequence_cleanup(bg_canonical_transcript, aa_mutation)
			if aa_change and aa_change != "":
				cons.append("%s:%s"%(consequence, aa_change))
			else:
				cons.append(consequence)

	return [impact_estimate, ";".join(cons)]


#########################################
def protein_coding_genes(cursor):
	standard_chromosomes = [str(i) for i in range(23)] + ['X','Y']
	genes = []
	chrom = {}
	qry  = "select approved_symbol, chromosome from icgc.hgnc "
	qry += "where locus_group='protein-coding gene'"
	for gene,chr in hard_landing_search(cursor,qry):
		if not chr in standard_chromosomes: continue
		genes.append(gene)
		chrom[gene] = chr

	return genes, chrom

#########################################
def co_ocurrence_raw(cursor, somatic_table, gene1, gene2):

	qry =  "select g1.gene_symbol,  g2.gene_symbol, s1.icgc_donor_id, s1.submitted_sample_id  "
	qry += "from mutation2gene g1, mutation2gene g2,  %s s1,  %s s2  " % (somatic_table, somatic_table)
	qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
	qry += "and s1.icgc_mutation_id=g1.icgc_mutation_id and g1.gene_symbol='%s' " % gene1
	qry += "and s2.icgc_mutation_id=g2.icgc_mutation_id and g2.gene_symbol='%s' " % gene2
	qry += "and s1.pathogenicity_estimate=1 and s1.reliability_estimate=1  "
	qry += "and s2.pathogenicity_estimate=1 and s2.reliability_estimate=1 "
	return search_db(cursor,qry)

#########################################
def quotify(something):
	if not something:
		return ""
	return "'{}'".format(something)

#########################################
# beats me, but using the intermediate (mutation2gene) is faster
# using inner join or not makes not difference
def co_ocurrence_w_group_count_inner_join(cursor, somatic_table, gene1, other_genes):
	qry =  "select count(distinct s1.icgc_donor_id) ct "
	qry += "from  %s s1,  %s s2  " % (somatic_table, somatic_table)
	#qry += "from  %s s1  " % (somatic_table)
	qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
	qry += "and s1.gene_symbol='%s' " % gene1
	#qry += "inner join %s as s2 on s1.icgc_donor_id=s2.icgc_donor_id " %somatic_table
	#qry += "where  s1.gene_symbol='%s' " % gene1
	group_string = (",".join([quotify(gene2) for gene2 in other_genes]))
	qry += "and s2.gene_symbol in (%s) " % group_string
	qry += "and s1.pathogenicity_estimate=1 and s1.reliability_estimate=1 "
	qry += "and s2.pathogenicity_estimate=1 and s2.reliability_estimate=1 "

	ret = search_db(cursor,qry)

	if not ret:
		search_db(cursor,qry,verbose=True)
		exit()
	return ret[0][0]

#########################################
def co_ocurrence_w_group_count(cursor, somatic_table, gene1, other_genes):
	qry =  "select count(distinct s1.icgc_donor_id) ct "
	qry += "from mutation2gene g1, mutation2gene g2,  %s s1,  %s s2  " % (somatic_table, somatic_table)
	qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
	qry += "and s1.icgc_mutation_id=g1.icgc_mutation_id and g1.gene_symbol='%s' " % gene1
	group_string = (",".join([quotify(gene2) for gene2 in other_genes]))
	qry += "and s2.icgc_mutation_id=g2.icgc_mutation_id and g2.gene_symbol in (%s) " % group_string
	qry += "and s1.pathogenicity_estimate=1 and s1.reliability_estimate=1 "
	qry += "and s2.pathogenicity_estimate=1 and s2.reliability_estimate=1 "

	ret = search_db(cursor,qry)

	if not ret:
		search_db(cursor,qry,verbose=True)
		exit()
	return ret[0][0]

#########################################
def create_gene_view(cursor, somatic_table, gene):
	view_name = "view_{}_{}_{}".format(somatic_table, gene, os.getpid())
	qry = "create view %s " % view_name
	qry += "as select distinct icgc_donor_id from  %s  " % somatic_table
	qry += "where gene_symbol='%s' " % gene
	qry += "and pathogenicity_estimate=1 and reliability_estimate=1"
	error_intolerant_search(cursor, qry)
	return view_name

#########################################
def create_gene_temp(cursor, somatic_table, gene):
	temp_name = "temp_{}_{}_{}".format(somatic_table, gene, os.getpid())
	qry = "create temporary table %s " % temp_name
	qry += "as select distinct icgc_donor_id from  %s  " % somatic_table
	qry += "where gene_symbol='%s' " % gene
	qry += "and pathogenicity_estimate=1 and reliability_estimate=1"
	error_intolerant_search(cursor, qry)
	return temp_name
#########################################
def create_gene_faux_temp(cursor, somatic_table, gene):
	temp_name = "temp_{}_{}_{}".format(somatic_table, gene, os.getpid())
	qry = "create  table %s " % temp_name
	qry += "as select distinct icgc_donor_id from  %s  " % somatic_table
	qry += "where gene_symbol='%s' " % gene
	qry += "and pathogenicity_estimate=1 and reliability_estimate=1"
	error_intolerant_search(cursor, qry)
	return temp_name
############
def drop_view(cursor, view_name):
	qry = "drop view if exists %s " % view_name
	error_intolerant_search(cursor, qry)
	return

# ./icgc_utils/kernprof.py -l  ....py
# python3 -m line_profiler .....py.lprof
# @profile
#########################################
def co_occurrence_count_with_a_temp(cursor, somatic_table, gene1_view, gene2):

	qry =  "select count(distinct s1.icgc_donor_id) ct "
	qry += "from  %s s1,  %s s2  " % (gene1_view, somatic_table)
	qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
	qry += "and s2.gene_symbol='%s' " % gene2
	qry += "and s2.pathogenicity_estimate=1 and s2.reliability_estimate=1 "
	ret = search_db(cursor,qry,verbose=True)
	if not ret:
		search_db(cursor,qry,verbose=True)
		retval = 0
	else:
		retval = ret[0][0]

	return retval

def co_occurrence_count(cursor, somatic_table, gene1, gene2):

	qry =  "select count(distinct s1.icgc_donor_id) ct "
	qry += "from  %s s1,  %s s2  " % (somatic_table, somatic_table)
	qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
	qry += "and s1.gene_symbol='%s' " % gene1
	qry += "and s2.gene_symbol='%s' " % gene2
	qry += "and s1.pathogenicity_estimate=1 and s1.reliability_estimate=1 "
	qry += "and s2.pathogenicity_estimate=1 and s2.reliability_estimate=1 "

	ret = search_db(cursor,qry,verbose=True)
	if not ret:
		search_db(cursor,qry,verbose=True)
		retval = 0
	else:
		retval = ret[0][0]

	return retval


#########################################
def patients_per_gene_breakdown(cursor, table):

	# this hinges on s.icgc_mutation_id=g.icgc_mutation_id
	# having icgc_mutation_id indexed both on s and g:
	qry  = "select g.gene_symbol symbol, count(distinct  s.icgc_donor_id) ct "
	qry += "from mutation2gene g, %s s  " % table
	qry += "where s.icgc_mutation_id=g.icgc_mutation_id and s.pathogenicity_estimate=1  "
	qry += "and s.reliability_estimate=1 "
	qry += "group by symbol"
	ret = hard_landing_search(cursor,qry)
	return dict(ret)

########################################
def genes_per_patient_breakdown(cursor, table):
	qry  = "select distinct s.icgc_donor_id donor, count(distinct g.gene_symbol) ct "
	qry += "from mutation2gene g, %s s  " % table
	qry += "where s.icgc_mutation_id=g.icgc_mutation_id and s.pathogenicity_estimate=1  "
	qry += "and s.reliability_estimate=1 "
	qry += "group by donor"
	ret = hard_landing_search(cursor,qry)
	return dict(ret)


#########################################
def patients_with_muts_in_gene_group(cursor, table, gene_list):

	# this hinges on s.icgc_mutation_id=g.icgc_mutation_id
	# having icgc_mutation_id indexed both on s and g:
	qry  = "select count(distinct  s.icgc_donor_id) ct "
	qry += "from mutation2gene g, %s s  " % table
	qry += "where s.icgc_mutation_id=g.icgc_mutation_id and s.pathogenicity_estimate=1  "
	qry += "and s.reliability_estimate=1 "
	group_string = (",".join([quotify(gene2) for gene2 in gene_list]))
	qry += "and g.gene_symbol in (%s)" %  group_string

	ret = hard_landing_search(cursor, qry)
	return ret[0][0]

########################################
def find_chromosome(cursor, gene):
	qry = "select chromosome from icgc.hgnc where approved_symbol = '%s'" % gene
	ret = search_db(cursor,qry)
	if not ret or ret ==[]:
		print("chromosome not found for %s (?)"%gene)
		search_db(cursor,qry,verbose=True)
		exit()
	return ret[0][0]

########################################
def get_donors(cursor, table):
	qry  = "select distinct(icgc_donor_id) from %s " % table
	return [ret[0] for ret in search_db(cursor,qry)]

def get_mutations(cursor, table, chromosome=None):
	qry = "select  distinct(icgc_mutation_id)  from %s " % table
	if chromosome: qry += "where chromosome='%s' " % chromosome
	ret = search_db(cursor,qry)
	if not ret: return []
	return [r[0] for r in ret]

def get_number_of_path_mutations_per_specimen(cursor, table, specimen_id):
	qry  = "select count(distinct icgc_mutation_id)  from %s " % table
	qry += "where  icgc_specimen_id = '%s' " % specimen_id
	qry += "and pathogenicity_estimate=1 and reliability_estimate=1 "
	return search_db(cursor,qry)[0][0]

def get_consequence(cursor, chromosome, mutation):
	qry  = "select consequence, aa_mutation from mutations_chrom_%s " % chromosome
	qry += "where icgc_mutation_id='%s' " % mutation
	return search_db(cursor,qry)[0]

def get_specimens_from_donor(cursor, table, icgc_donor_id):
	qry = "select  distinct(icgc_specimen_id)  from %s " % table
	qry += "where icgc_donor_id = '%s'" % icgc_donor_id
	return [r[0] for r in search_db(cursor,qry)]

def get_specimen_type(cursor, tumor_short, spec_ids):
	specimen_type = {}
	for spec_id in spec_ids:
		qry = " select specimen_type from %s_specimen " % tumor_short
		qry += "where icgc_specimen_id = '%s'" % spec_id
		specimen_type[spec_id] = search_db(cursor,qry)[0][0]
	return specimen_type

def get_mutations_from_donor(cursor, table, icgc_donor_id):
	qry = "select  distinct(icgc_mutation_id)  from %s " % table
	qry += "where icgc_donor_id = '%s'" % icgc_donor_id
	return [r[0] for r in search_db(cursor,qry)]

def mutation_provenance(cursor, table, icgc_donor_id, icgc_mutation_id):
	qry = "select  distinct(icgc_specimen_id)  from %s " % table
	qry += "where icgc_donor_id='%s' and icgc_mutation_id='%s'" % (icgc_donor_id, icgc_mutation_id)
	return [r[0] for r in search_db(cursor,qry)]

def mutation_count_per_donor(cursor, somatic_table):
	qry = "select icgc_donor_id, count(*) as c  from %s " % somatic_table
	qry += "where pathogenicity_estimate=1 and reliability_estimate=1 "
	qry += "group by  icgc_donor_id"
	ret = hard_landing_search(cursor, qry)
	return dict(ret)

def get_number_of_genes_affected(cursor, table):
	qry = "select  count(distinct gene_symbol) from mutation2gene where icgc_mutation_id in "
	qry += "(select distinct icgc_mutation_id from %s " % table
	qry += "where pathogenicity_estimate=1 and reliability_estimate=1)"
	return hard_landing_search(cursor,qry)[0][0]

#########################################
def mutations_in_gene_old(cursor, approved_symbol):
	qry  = "select ensembl_gene_id_by_hgnc, ensembl_gene_id, chromosome from hgnc "
	qry += "where approved_symbol = '%s'" % approved_symbol
	ensembl_gene_id_by_hgnc, ensembl_gene_id, chromosome = search_db(cursor,qry)[0]
	if ensembl_gene_id_by_hgnc != ensembl_gene_id:
		print("Ensembl id mismatch: (ensembl_gene_id_by_hgnc, ensembl_gene_id)")
		print(approved_symbol, ensembl_gene_id_by_hgnc, ensembl_gene_id)
		exit()

	qry  = "select m.icgc_mutation_id from mutations_chrom_%s m, locations_chrom_%s l "  % (chromosome, chromosome)
	qry += "where m.pathogenicity_estimate=1 and m.start_position=l.position "
	qry += "and l.gene_relative like '%%%s%%' " % ensembl_gene_id
	ret = search_db(cursor,qry, verbose=True)
	if not ret: return []
	return [r[0] for r in ret]

#########################################
def mutations_in_gene(cursor, approved_symbol):
	qry  = "select icgc_mutation_id from mutation2gene "
	qry += "where gene_symbol='%s'" % approved_symbol
	ret = search_db(cursor,qry, verbose=True)
	if not ret: return []
	return [r[0] for r in ret]

#########################################
def pathogenic_mutations_in_gene(cursor, approved_symbol, chromosome, use_reliability=True):
	qry  = "select map.icgc_mutation_id from mutation2gene map, mutations_chrom_%s mut " % chromosome
	qry += "where map.gene_symbol='%s' " % approved_symbol
	qry += "and map.icgc_mutation_id=mut.icgc_mutation_id "
	qry += "and mut.pathogenicity_estimate=1 "
	if use_reliability: qry += "and mut.reliability_estimate=1 "
	ret = search_db(cursor,qry, verbose=True)
	if not ret: return []
	return [r[0] for r in ret]


#########################################
def ensembl_gene_id2approved_symbol(cursor, ensembl_gene_id):
	symbol = None
	qry = "select approved_symbol from icgc.hgnc where ensembl_gene_id='%s'"% ensembl_gene_id
	ret = search_db(cursor,qry)
	if not ret:
		new_id =  attempt_resolve_deprecated(cursor, ensembl_gene_id)
		if new_id:
			qry = "select approved_symbol from icgc.hgnc where ensembl_gene_id='%s'"% new_id
			ret = search_db(cursor,qry)
			if ret: symbol = ret[0][0]

	else:
		symbol = ret[0][0]
	# if not resolved return the original ensebl_gene_id
	return symbol if symbol else ensembl_gene_id


#########################################
def attempt_resolve_deprecated(cursor, stable_id_old, verbose=False):
	qry = "select new_id from icgc.ensembl_deprecated_ids where old_id='%s'" % stable_id_old
	ret = search_db(cursor,qry,verbose=verbose)
	if not ret: return None
	return ret[0][0]

# ./icgc_utils/kernprof.py -l 39_reannotate_missense_mutations.py
# python3 -m line_profiler 39_reannotate_missense_mutations.py.lprof
#@profile
#########################################
def gene_stable_id_2_canonical_transcript_id(cursor, gene_stable_id, verbose=False):
	qry  = "select  distinct(canonical_transcript) from icgc.ensembl_ids where  gene ='%s' " % gene_stable_id
	ret = search_db(cursor,qry)
	if not ret:
		new_id =  attempt_resolve_deprecated(cursor, gene_stable_id, verbose)
		if not new_id:
			if verbose: print("No canonical transcript and no new id found for %s "% gene_stable_id)
			return None
		qry  = "select  distinct(canonical_transcript) from icgc.ensembl_ids where  gene ='%s' " % new_id
		ret = search_db(cursor,qry)
		if ret:
			return ret[0][0]
		if verbose: print("No canonical transcript id found for %s, mapped to %s" % (gene_stable_id,new_id))
		return None
	elif len(ret) != 1:
		if verbose: print("No unique canonical transcript id found for %s" % gene_stable_id)
		return None
	return ret[0][0]

#########################################
def list_of_transcript_ids_2_canonical_transcript_id(cursor, list_of_stable_transcript_ids):
	# list_od_stable_transcript_ids - refers to ensembl --> ENST00... identifier
	ensts = ",".join(["'%s'"%enst for enst in list_of_stable_transcript_ids])
	qry  = "select distinct(canonical_transcript) from ensembl_ids  where transcript in  (%s) " % ensts
	ret = search_db(cursor,qry)
	if not ret:
		#print("Warning: no canonical transcript found for %s" % ensts)
		#print("Qry was: ", qry)
		return []
	# there may be multiple canonical transcripts if the list of transcripts belongs to different genes
	return [r[0] for r in ret]


#########################################
def approved_symbol2ensembl_canonical_transcript(cursor, gene_symbol):
	qry  = "select ensembl_gene_id from icgc.hgnc "
	qry += "where approved_symbol = '%s' " % gene_symbol
	[ensembl_gene_id] = hard_landing_search(cursor,qry)[0]
	return gene_stable_id_2_canonical_transcript_id(cursor, ensembl_gene_id)

####################################################
def nonzerolen(stringlist):
	return list(filter(lambda x: len(x)>0, stringlist))
####################################################
def gene_group_cdna_length(cursor, gene_group):

	if len(gene_group)==0: return 0

	gene_string = ",".join([quotify(g) for g in gene_group])
	qry  = "select distinct(ensembl_gene_id) from hgnc "
	qry += "where approved_symbol in (%s)" % gene_string
	ensembl_gene_ids = nonzerolen([r[0] for r in hard_landing_search(cursor,qry)])

	ensid_string = ",".join([quotify(ensid) for ensid in ensembl_gene_ids])
	qry  = "select distinct(canonical_transcript) from ensembl_ids "
	qry += "where gene in (%s)" % ensid_string
	ensembl_transcript_ids = nonzerolen([r[0] for r in hard_landing_search(cursor,qry)])

	enstr_string = ",".join([quotify(ensid) for ensid in ensembl_transcript_ids])
	qry  = "select length(sequence) from ensembl_coding_seqs "
	qry += "where transcript_id in (%s)" % enstr_string
	length = sum([r[0] for r in hard_landing_search(cursor,qry)])

	return length
