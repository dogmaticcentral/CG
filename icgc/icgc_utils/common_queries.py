
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
from icgc_utils.mysql   import  *
from icgc_utils.CrossMap import *

#########################################
def gnomad_mutations (cursor, gene_symbol):

	mutations = []

	chromosome = find_chromosome(cursor, gene_symbol)
	#column_names
	colnames = get_column_names(cursor,"gnomad","gnomad_freqs_chr_1")

	# brute force approach seems to be fast enough for a single gene
	qry = "select * from gnomad.gnomad_freqs_chr_{} where consequences like '%|{}|%' ".format(chromosome, gene_symbol)
	ret = search_db(cursor,qry)
	if not ret:
		print("nothing found for {}, chromosome {}".format(gene_symbol, chromosome))
		exit()
	for line in ret:
		named_fields = dict(list(zip(colnames,line)))
		relevant_variants = []
		for description in named_fields['consequences'].split(","):
			if not 'RPL5' in description: continue
			if not 'missense' in description: continue
			description_field = description.split("|")
			# I don't have ensembl info here - in a more through implementation one should
			# at least go for annotator here
			# for now just hope that the uniprot is canonical
			if len(description_field[2])==0: continue
			relevant_variants.append(description)
		if len(relevant_variants)==0: continue
		if float(named_fields['variant_count'])<2: continue
		freqency = float(named_fields['variant_count'])/named_fields['total_count']
		#print "%.1e" % freqency,
		for description in relevant_variants:
			description_field = description.split("|")
			#print "  ", description_field[7], description_field[8], # example:  280 V/A
			mutations.append(description_field[8].split("/")[0] + description_field[7] + description_field[8].split("/")[1])
		#print

	return list(set(mutations))

#########################################
def translate_positions(positions, chromosome, from_assembly, to_assembly, rootname=None):

	if from_assembly == to_assembly:
		return positions
	# GRCh37 and hg19 only differ for MT
	if (from_assembly.lower() in ['grch37', 'hg19']) and (to_assembly.lower() in ['grch37', 'hg19']) and (chromosome != "MT"):
		return positions

	if not rootname: rootname="{}.{}.{}.{}".format(chromosome, from_assembly, to_assembly, str(os.getpid()))
	# otherwise we'll need  tools to translate
	chain_file ="/storage/databases/liftover/{}To{}.over.chain".format(from_assembly, to_assembly.capitalize())
	if not os.path.exists(chain_file):
		print(chain_file, "not found")
		exit()

	outfile = "%s.tsv"%rootname
	outf = open (outfile,"w")
	for p in positions:
		chr = chromosome if 'chr' in chromosome else 'chr'+chromosome
		outf.write("\t".join([chr, str(p), str(p)]) + "\n")
	outf.close()

	# this is CrossMap now
	outfile_translated  =  "%s.translated.tsv"%rootname
	(map_tree, target_chrom_sizes, source_chrom_sizes) = read_chain_file(chain_file, print_table=False)
	crossmap_bed_file(map_tree, outfile, outfile_translated)

	#read binding regions back in
	with open(outfile_translated,"r") as inf:
		new_positions = [line.split("\t")[1] for line in inf.read().split("\n") if len(line.replace(" ",""))>0]

	if len(new_positions) != len(positions):
		print("translation mismatch, {}, chromosome {}, {} to {}".format(os.getcwd(), chromosome, from_assembly, to_assembly))
		print("process", os.getpid(), "exiting")
		exit(1)
	else:
		# remove aux files
		os.remove(outfile)
		os.remove(outfile_translated)
		os.remove(outfile_translated+".unmap") # this file should probably checked - it should be empty

	return new_positions

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
		qry += "from %s_specimen where icgc_specimen_id in (%s)" % (
		tumor, ",".join(["'%s'" % id for id in entries_per_specimen.keys()]))
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
def aa_change_cleanup(cursor, aa_change):
	if not aa_change: return ""
	if aa_change=="": return aa_change
	change = {};
	for enst_change in aa_change.split(";"):
		[e, c] = enst_change.split(":")
		change[e] = c
	enst_canonical = list_of_transcript_ids_2_canonical_transcript_id(cursor, list(change.keys()))
	if not enst_canonical or not enst_canonical in change: return aa_change
	return change[enst_canonical]

#########################################
def find_53_status(cursor, tumor_short, specimen):
	# g = gene
	# m = mutation
	# v = variant
	qry  = "select g.icgc_mutation_id, v.pathogenicity_estimate, m.consequence, m.aa_mutation, l.transcript_relative "
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
		if line[-1]!=None and 'splice' in line[-1]:
			cons.append(line[-1])
			continue
		if line[2] and line[3]:
			aa_change = aa_change_cleanup(cursor, line[3])
			if aa_change and aa_change != "":
				cons.append("%s:%s"%(line[2], aa_change))
			else:
				cons.append(line[2])
	return [impact_estimate, ";".join(cons)]


#########################################
def protein_coding_genes(cursor):
	standard_chromosomes = [str(i) for i in range(23)] + ['X','Y']
	genes = []
	chrom = {}
	qry  = "select approved_symbol, chromosome from icgc.hgnc "
	qry += "where locus_group='protein-coding gene'"
	for gene,chr in search_db(cursor,qry):
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
	qry += "and s1.pathogenic_estimate=1 and s1.reliability_estimate=1  "
	qry += "and s2.pathogenic_estimate=1 and s2.reliability_estimate=1 "
	return search_db(cursor,qry)

#########################################
def quotify(something):
	if not something:
		return ""
	if type(something)==str:
		return "\'"+something+"\'"
	else:
		return str(something)


#########################################
def co_ocurrence_w_group_count(cursor, somatic_table, gene1, other_genes):
	qry =  "select count(distinct s1.icgc_donor_id) ct "
	qry += "from mutation2gene g1, mutation2gene g2,  %s s1,  %s s2  " % (somatic_table, somatic_table)
	qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
	qry += "and s1.icgc_mutation_id=g1.icgc_mutation_id and g1.gene_symbol='%s' " % gene1
	group_string = (",".join([quotify(gene2) for gene2 in other_genes]))
	qry += "and s2.icgc_mutation_id=g2.icgc_mutation_id and g2.gene_symbol in (%s) " % group_string
	qry += "and s1.pathogenic_estimate=1 and s1.reliability_estimate=1 "
	qry += "and s2.pathogenic_estimate=1 and s2.reliability_estimate=1 "

	ret = search_db(cursor,qry)

	if not ret:
		search_db(cursor,qry,verbose=True)
		exit()
	return ret[0][0]


#########################################
def co_ocurrence_count(cursor, somatic_table, gene1, gene2):

	if True: # this is still  twice as fast as the search below
		# are we running thruoght the same row twice bcs of s1 <-> s2?
		# still, distinct should get rid of double counting
		qry =  "select count(distinct s1.icgc_donor_id) ct "
		qry += "from mutation2gene g1, mutation2gene g2,  %s s1,  %s s2  " % (somatic_table, somatic_table)
		qry += "where s1.icgc_donor_id=s2.icgc_donor_id "
		qry += "and s1.icgc_mutation_id=g1.icgc_mutation_id and g1.gene_symbol='%s' " % gene1
		qry += "and s2.icgc_mutation_id=g2.icgc_mutation_id and g2.gene_symbol='%s' " % gene2
		qry += "and s1.pathogenic_estimate=1 and s1.reliability_estimate=1  "
		qry += "and s2.pathogenic_estimate=1 and s2.reliability_estimate=1 "
		ret = search_db(cursor,qry)
		if not ret:
			search_db(cursor,qry,verbose=True)
			exit()
		return ret[0][0]
	else:
		donors = {}
		for gene in [gene1, gene2]:
			qry  = "select distinct s.icgc_donor_id "
			qry += "from mutation2gene g,  %s s  " % (somatic_table)
			qry += "where s.icgc_mutation_id=g.icgc_mutation_id and g.gene_symbol='%s' " % gene
			qry += "and s.pathogenic_estimate=1 and s.reliability_estimate=1  "
			ret  = search_db(cursor,qry)
			if not ret:
				search_db(cursor,qry,verbose=True)
				exit()
			donors[gene] = [r[0] for r in ret]

		return len(set(donors[gene1])&set(donors[gene2]))

#########################################
def patients_per_gene_breakdown(cursor, table):

	# this hinges on s.icgc_mutation_id=g.icgc_mutation_id
	# having icgc_mutation_id indexed both on s and g:
	qry  = "select g.gene_symbol symbol, count(distinct  s.icgc_donor_id) ct "
	qry += "from mutation2gene g, %s s  " % table
	qry += "where s.icgc_mutation_id=g.icgc_mutation_id and s.pathogenic_estimate=1  "
	qry += "and s.reliability_estimate=1 "
	qry += "group by symbol"
	ret = search_db(cursor,qry)
	if not ret:
		search_db(cursor,qry, verbose=True)
		exit()
	return dict(ret)

#########################################
def patients_with_muts_in_gene_group(cursor, table, gene_list):

	# this hinges on s.icgc_mutation_id=g.icgc_mutation_id
	# having icgc_mutation_id indexed both on s and g:
	qry  = "select count(distinct  s.icgc_donor_id) ct "
	qry += "from mutation2gene g, %s s  " % table
	qry += "where s.icgc_mutation_id=g.icgc_mutation_id and s.pathogenic_estimate=1  "
	qry += "and s.reliability_estimate=1 "
	group_string = (",".join([quotify(gene2) for gene2 in gene_list]))
	qry += "and g.gene_symbol in (%s)" %  group_string

	ret = search_db(cursor,qry)
	if not ret:
		search_db(cursor,qry, verbose=True)
		exit()
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
def try_to_resolve(cursor, old_ensembl_gene_id):
	# not sure how stable or reliable this is, but if there is no common transcript,
	# something iss seriously foul with the old_ensembl_gene_id
	switch_to_db(cursor,"homo_sapiens_core_91_38")
	qry = "select distinct(gene_stable_id) from gene_archive "
	qry += "where transcript_stable_id in  "
	qry +="(select transcript_stable_id from gene_archive where gene_stable_id = '%s')" % old_ensembl_gene_id
	ret = search_db(cursor,qry)
	if not ret: return None

	candidate_ids = [r[0] for r in ret if r[0]!=old_ensembl_gene_id]
	if len(candidate_ids)==0: return None  # there should be another identifier, besides the one we started from

	latest_ensembl_entries = []
	for  candidate in candidate_ids:
		qry = "select * from gene  where stable_id = '%s'" % candidate
		if not search_db(cursor,qry): continue  # this is another old identifier
		latest_ensembl_entries.append(candidate)
	#I'm not sure what to make of this if there are two live identifiers
	# that the old one maps to
	if len(latest_ensembl_entries) != 1: return None

	new_ensembl_gene_id = latest_ensembl_entries[0]

	qry = "select approved_symbol from icgc.hgnc where ensembl_gene_id='%s'"% new_ensembl_gene_id
	ret = search_db(cursor,qry)
	if not ret: return None
	return ret[0][0]


def get_approved_symbol(cursor, ensembl_gene_id):
	qry = "select approved_symbol from hgnc where ensembl_gene_id='%s'"% ensembl_gene_id
	ret = search_db(cursor,qry)
	if not ret:
		symbol = try_to_resolve(cursor, ensembl_gene_id)
		# if it cannot be resolved, just use the ensembl_id_itself
		if not symbol: symbol=ensembl_gene_id
		switch_to_db(cursor,"icgc")
	else:
		symbol = ret[0][0]
	return symbol

#########################################
def gene_stable_id_2_canonical_transcript_id(cursor, gene_stable_id, verbose=False):
	qry  = "select  distinct(canonical_transcript) from icgc.ensembl_ids where  gene ='%s' " % gene_stable_id
	ret = search_db(cursor,qry)
	if not ret or len(ret) != 1:
		if verbose: print("Warning: no unique canonical id could be found for %s" % gene_stable_id)
		return None
	return ret[0][0]


def list_of_transcript_ids_2_canonical_transcript_id(cursor, list_of_stable_transcript_ids):
	# list_od_stable_transcript_ids - refers to ensembl --> ENST00... identifier
	ensts = ",".join(["'%s'"%enst for enst in list_of_stable_transcript_ids])
	qry  = "select distinct(canonical_transcript) from ensembl_ids  where transcript in  (%s) " % ensts
	ret = search_db(cursor,qry)
	if not ret or len(ret) != 1:
		print("Warning: no unique canonical transcript id could be found for %s" % ensts)
		print("Qry was: ", qry)
		return None
	return ret[0][0]


