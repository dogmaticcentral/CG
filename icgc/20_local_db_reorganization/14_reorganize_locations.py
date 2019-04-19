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

import time

from config import Config
from icgc_utils.processes import *
from icgc_utils.annovar import *
from icgc_utils.CrossMap import *


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
	outfname = "%s.%s.avinput" % (variants_table,assembly)
	outf = open(outfname, 'w')
	for position in positions:
		outrow = "\t".join([chromosome, str(position), str(position), '-', 'A'])
		outf.write(outrow+"\n")
	outf.close()
	return outfname


#########################################
def translate_positions(positions, chromosome, from_assembly, to_assembly, rootname=None):

	if from_assembly == to_assembly:
		return positions
	# GRCh37 and hg19 only differ for MT
	if (from_assembly.lower() in ['grch37', 'hg19']) and (to_assembly.lower() in ['grch37', 'hg19']) and (chromosome != "MT"):
		return positions

	if not rootname: rootname=str(os.getpid())
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
def make_location_tsv(cursor, avoutput):

	inf  = open (avoutput, "r")
	header_fields = None
	tsv_list = []
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
		tsv_list.append("\t".join([str(position), gene_relative, transcript_relative]))
	inf.close()
	return tsv_list

# ./icgc_utils/kernprof.py -l 14_reorganize_locations.py
# python3 -m line_profiler 14_reorganize_locations.py.lprof
# @profile
#########################################
def reorganize_locations(cursor, ref_assembly, chromosome, tables):

	time0 = time.time()
	avinputs = []
	for variants_table  in tables:
		# which assemblies do we have in this story (as of Apr 2019 it was only GRCh37 - takes 6 mins to find this out)
		assemblies = ['GRCh37']
		# qry = "select distinct assembly  from %s" % variants_table
		# assemblies = [r[0] for r in search_db (cursor, qry)]
		for assembly in assemblies:
			positions = get_positions(cursor, variants_table, chromosome, assembly)
			if len(positions)==0: continue
			positions_translated = translate_positions(positions, chromosome, assembly, ref_assembly)
			# write_annovar_input will fill avinput dictionary
			avinput  = write_annovar_input(positions_translated, variants_table, ref_assembly, chromosome)
			avinputs.append(avinput)
			#print("\t\t\t", variants_table, chromosome, assembly, "annovar input written")
	if len(avinputs)==0:
		print("pid %d: no annovar input written (?)" % os.getpid())
		return
	# concatenate all input, remove duplicates
	subprocess.call(["bash","-c", " cat %s | sort -gk2 | uniq > all.avinput" % (" ".join(avinputs)) ])
	# the positions in the avinput are already hg19
	avoutput = run_annovar("all.avinput", ref_assembly, "all", chromosome=="MT")
	time1 = time.time()
	print("\t\t\t chromosome %s avoutput written to %s in %.3f mins "%(chromosome, avoutput, float(time1-time0)/60))

	tsv_list = make_location_tsv(cursor, avoutput)
	time2 = time.time()
	print("\t\t\t chromosome %s tsv list created in %.3f mins" % (chromosome, float(time2-time1)/60))

	tsv_name = "locations_chrom_%s.tsv" % chromosome
	outf = open (tsv_name,"w")
	for entry in tsv_list:
		outf.write(entry + "\n")
	outf.close()
	time3 = time.time()
	print("\t\t\t chromosome %s tsv written to %s in %.3f mins" % (chromosome, tsv_name, float(time3-time2)/60))

	# eat this
	qry = "load data local infile '%s' into table locations_chrom_%s" % (tsv_name, chromosome)
	search_db(cursor,qry,verbose=True)
	time4 = time.time()
	print("\t\t\t chromosome %s loaded in %.3f mins"%(chromosome, float(time4-time3)/60))

	return


#########################################
def reorganize(chromosomes, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	ref_assembly = other_args[0]
	tables = other_args[1]
	home = os.getcwd()

	for chrom in chromosomes:

		# make a workdir and move there
		workdir  = "chrom_%s" % chrom
		workpath = "{}/annovar/{}".format(home,workdir)
		if not os.path.exists(workpath): os.makedirs(workpath)
		os.chdir(workpath)

		time0 = time.time()
		print("====================")
		print("reorganizing locations for chrom %s, pid %d" % (chrom, os.getpid()))
		reorganize_locations(cursor, ref_assembly, chrom, tables)
		time1 = time.time()
		print(("chrom %s  done in %.3f mins" % (chrom, float(time1-time0)/60)), os.getpid())

	cursor.close()
	db.close()

	return

#########################################
#########################################
#########################################
def main():

	#print("Disabled.")
	#exit()

	ref_assembly = 'hg19' # this is the assembly I would like to see for all coords in location tables
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	chromosomes = [str(i) for i in range(1,23)] + ["X", "Y"]
	number_of_chunks = 12  # myISAM does not deadlock

	#chromosomes = ["Y"]
	#number_of_chunks = 1
	parallelize(number_of_chunks, reorganize, chromosomes, [ref_assembly, tables], round_robin=True)


#########################################
if __name__ == '__main__':
	main()
