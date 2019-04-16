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
import time

from config import Config
from icgc_utils.common_queries import  *
from icgc_utils.processes import  *
from icgc_utils.annovar import *
from icgc_utils.CrossMap import *
from subprocess import PIPE

variant_columns = ['icgc_mutation_id', 'chromosome','icgc_donor_id', 'icgc_specimen_id', 'icgc_sample_id',
					'submitted_sample_id','control_genotype', 'tumor_genotype', 'total_read_count', 'mutant_allele_read_count']

#  'aa_mutation',  'consequence_type', and 'pathogenic_estimate'  will be filled separately
mutation_columns = ['icgc_mutation_id', 'start_position', 'end_position', 'assembly',
					'mutation_type', 'mutated_from_allele', 'mutated_to_allele', 'reference_genome_allele']

location_columns = ['position', 'gene_relative', 'transcript_relative']

################################################################
# stop_retained: A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
consequence_vocab = ['stop_lost', 'synonymous', 'inframe_deletion', 'inframe_insertion', 'stop_gained',
					 '5_prime_UTR_premature_start_codon_gain',
					 'start_lost', 'frameshift', 'disruptive_inframe_deletion', 'stop_retained',
					 'exon_loss', 'disruptive_inframe_insertion', 'missense']

# location_vocab[1:4] is gene-relative
# location_vocab[4:] is transcript-relative
location_vocab = ['intergenic_region', 'intragenic', 'upstream', 'downstream',
				  'intron', 'exon',  'coding_sequence', '5_prime_UTR',  '3_prime_UTR',
				  'initiator_codon', 'splice_donor','splice_acceptor', 'splice_region'
				  ]

# this is set literal
pathogenic = {'missense', 'exon_loss',
				'stop_lost', 'stop_gained', 'start_lost',
				'5_prime_UTR_premature_start_codon_gain',
				'frameshift', 'inframe_deletion', 'inframe_insertion',
				'disruptive_inframe_deletion', 'disruptive_inframe_insertion',
				'inframe' }  # there is no way we can know at this level whether an inframe change is nondsisruptive
									# more likely it is than not

consequence2location = {
	'synonymous':'exon',
	'stop_retained':'exon',
	'missense':'exon',
	'exon_loss':'exon',
	'stop_lost':'exon',
	'stop_gained':'exon',
	'start_lost':'exon',
	'5_prime_UTR_premature_start_codon_gain':'5_prime_UTR',
	'frameshift':'exon',
	'inframe_deletion':'exon',
	'inframe_insertion':'exon',
	'disruptive_inframe_deletion':'exon',
	'disruptive_inframe_insertion':'exon',
	'inframe':'exon'
}

benign = {'synonymous', 'stop_retained'}


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
def get_positions(cursor, variants_table, chromosome, assembly):
	positions = set()
	for position_significance in ['start', 'end']:
		qry =  "select distinct %s_position from %s  " % (position_significance, variants_table)
		qry += "where chromosome='%s' and assembly='%s' " %(chromosome, assembly)
		qry += "and gene_affected is not null and gene_affected !='' "
		ret  = search_db (cursor, qry)
		if ret:  positions |= set([r[0] for r in ret])
	return positions


#########################################
def write_annovar_input(positions, variants_table,assembly,chromosome):
	# we are duping annovar into giving us the location relative to each transcript
	# we are not interested in annotation (not here) therefore we give
	# the annovar some generic variant - like insert 'A' at given position
	# annovar input file(s)
	outfname = "%s.%s.%s.avinput" % (variants_table,assembly,chromosome)
	outf = open(outfname, 'w')
	for position in positions:
		outrow = "\t".join([chromosome, str(position), str(position), '-', 'A'])
		outf.write(outrow+"\n")
	outf.close()
	subprocess.call(["bash","-c", "sort %s | uniq > %s.tmp" % (outfname, outfname)])
	os.rename(outfname+".tmp", outfname)
	return outfname


#########################################
def translate_positions(positions, chromosome, from_assembly, to_assembly, rootname):

	if from_assembly == to_assembly:
		return positions

	# otherwise we'll need  tools to translate
	chain_file ="/storage/databases/liftover/{}To{}.over.chain".format(from_assembly, to_assembly.capitalize())
	if not os.path.exists(chain_file):
		print(chain_file, "not found")
		exit()

	outfile = "%s.tsv"%rootname
	outf = open (outfile,"w")
	for p in positions:
		outf.write("\t".join([chromosome, str(p), str(p)]) + "\n")
	outf.close()

	# this is CrossMap now
	outfile_translated  =  "%s.translated.tsv"%rootname
	(map_tree, target_chrom_sizes, source_chrom_sizes) = read_chain_file(chain_file, print_table = False)
	crossmap_bed_file(map_tree, outfile, outfile_translated)

	#read binding regions back in
	with open(outfile_translated,"r") as inf:
		new_positions = [line.split("\t")[1] for line in inf.read().split("\n") if len(line.replace(" ",""))>0]
	# remove aux files
	os.remove(outfile)
	os.remove(outfile_translated)
	os.remove(outfile_translated+".unmap") # this file should probably checked - it should be empty

	return new_positions


#########################################
def store_location_info(cursor, chromosome, avoutput):

	location_table = "locations_chrom_%s" % chromosome
	location_columns = ['position', 'gene_relative', 'transcript_relative']
	inf = open (avoutput, "r")
	header_fields = None
	for line in inf:
		if not header_fields:
			header_fields = parse_annovar_header_fields(line)
			if not header_fields:
				print("Error parsing annovar: no header line found")
				exit()
			continue
		[start, end, gene_relative, transcript_relative] = parse_annovar_line(cursor, header_fields,  line)
		position = start
		if start!=end:
			print("why is start not == end here?")
			exit()
		location_values = [str(position), quotify(gene_relative), quotify(transcript_relative)]
		#print(location_table, location_values)
		insert(cursor, location_table, location_columns, location_values)


#  ./icgc_utils/kernprof.py -l 14_reorganize_locations.py
# python3 -m line_profiler 14_reorganize_locations.py.lprof
# @profile
#########################################
def reorganize_locations(cursor, ref_assembly, variants_table):

	# which assemblies do we have in this story
	qry = "select distinct assembly  from %s" % variants_table
	assemblies = [r[0] for r in search_db (cursor, qry)]
	print(variants_table, "assemblies:",assemblies)

	chromosomes = [str(i) for i in range(1,23)] + ["X", "Y"]
	for chromosome in chromosomes:
		print("\t", variants_table, chromosome)
		rootname = "%s.%s" % (variants_table,chromosome)
		for assembly in assemblies:
			positions = get_positions(cursor, variants_table, chromosome, assembly)
			if len(positions)==0: continue
			positions_translated = translate_positions(positions, chromosome, assembly, ref_assembly, rootname)
			# write_annovar_input will fill avinput dictionary
			avinput  = write_annovar_input(positions_translated, variants_table, ref_assembly, chromosome)
			# the positions in the avinput are already hg19
			avoutput = run_annovar(avinput, ref_assembly, rootname, chromosome=="MT")
			store_location_info(cursor, chromosome, avoutput)
			#print(avoutput)


#########################################
def reorganize(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	ref_assembly = other_args[0]

	home = os.getcwd()

	for table in tables:

		# make a workdir and move there
		workdir  = table + "_annovar"
		workpath = "{}/annovar/{}".format(home,workdir)
		if not os.path.exists(workpath): os.makedirs(workpath)
		os.chdir(workpath)

		time0 = time.time()
		print("====================")
		print("reorganizing locations from", table, os.getpid())
		reorganize_locations(cursor, ref_assembly, table)
		time1 = time.time()
		print(("\t\t %s  done in %.3f mins" % (table, float(time1-time0)/60)), os.getpid())

	cursor.close()
	db.close()

	return

#########################################
#########################################
#########################################
def main():
	ref_assembly = 'hg19' # this is the assembly I would like to see all coords in location tables
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	number_of_chunks = 12  # myISAM does not deadlock
	parallelize(number_of_chunks, reorganize, tables, [ref_assembly])


#########################################
if __name__ == '__main__':
	main()
