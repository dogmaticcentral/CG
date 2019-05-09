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
from icgc_utils.utils import *

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
def write_positions(positions, variants_table,assembly,chromosome):
	# we are duping annovar into giving us the location relative to each transcript
	# we are not interested in annotation (not here) therefore we give
	# the annovar some generic variant - like insert 'A' at given position
	# annovar input file(s)
	outfname = "%s.txt" % (variants_table)
	outf = open(outfname, 'w')
	outf.write("\n".join([str(p) for p in positions])+"\n")
	outf.close()
	return outfname


#########################################
def describe_position(cursor, chrom, position, target_range):

	position -= 1 # UCSC labels positions from 0
	qry = "select ens_gene_id, ens_transcript_id, cds_start, cds_end, exon_count, exon_starts, exon_ends "
	qry += "from coords_chrom_%s " % chrom
	qry += "where tx_start between %d and %d " %(position-target_range, position+target_range)
	qry += "and %d between tx_start and tx_end " % position # 2 min

	ret = search_db (cursor, qry)
	if not ret:
		return [None, None]
	gene_relative_set = set()
	transcript_relative_list = []
	# the following expects that exon positions are sorted in the increasing order
	for line in ret:
		[ens_gene_id, ens_transcript_id, cds_start, cds_end, exon_count, exon_starts, exon_ends] = line
		gene_relative_set.add(ens_gene_id)
		if position<cds_start or position>cds_end:
			transcript_relative_list.append("{}:{}".format(ens_transcript_id, "UTR"))
		else:
			# ucscs exon list may end in ",' but not sure if it is quaranteed
			estart = [int(p) for p in exon_starts.split(",") if len(p)>0]
			eend = [int(p) for p in exon_ends.split(",") if len(p)>0]
			for i in range(exon_count):
				if position<estart[i]:
					if position>estart[i]-4:
						transcript_relative_list.append("{}:{}".format(ens_transcript_id, "splice"))
					else:
						transcript_relative_list.append("{}:{}".format(ens_transcript_id, "intronic"))
					break
				if position<=eend[i]+3:
					if position>eend[i]:
						transcript_relative_list.append("{}:{}".format(ens_transcript_id, "splice"))
					else:
						transcript_relative_list.append("{}:{}".format(ens_transcript_id, "exonic"))
					break
	return [";".join(list(gene_relative_set)), ";".join(transcript_relative_list)]


#########################################
def  longest_gene_range(cursor, chrom):
	qry = "select tx_end-tx_start as length from coords_chrom_%s order by length desc limit 1" % chrom
	return search_db(cursor, qry)[0][0]+1

#########################################
def make_location_tsv(cursor, chrom, pos_file):

	target_range = longest_gene_range(cursor, chrom)

	inf  = open (pos_file, "r")
	tsv_list = []
	for line in inf:
		position = int(line.rstrip())
		gene_relative, transcript_relative = describe_position(cursor, chrom, position, target_range)
		if gene_relative:
			tsv_list.append("\t".join([str(position), gene_relative, transcript_relative]))
	inf.close()
	return tsv_list

# ./icgc_utils/kernprof.py -l 14_reorganize_locations.py
# python3 -m line_profiler 14_reorganize_locations.py.lprof
# @profile
#########################################
def reorganize_locations(cursor, ref_assembly, chromosome, somatic_temp_tables):

	time0 = time.time()
	pos_files = []
	for variants_table  in somatic_temp_tables:
		# which assemblies do we have in this story (as of Apr 2019 it was only GRCh37 - takes 6 mins to find this out)
		assemblies = ['GRCh37']
		# qry = "select distinct assembly  from %s" % variants_table
		# assemblies = [r[0] for r in search_db (cursor, qry)]
		for assembly in assemblies:
			positions = get_positions(cursor, variants_table, chromosome, assembly)
			if len(positions)==0:continue
			positions_translated = translate_positions(positions, chromosome, assembly, ref_assembly)
			pos_file  = write_positions(positions_translated, variants_table, ref_assembly, chromosome)
			pos_files.append(pos_file)

			print("\t\t\t", variants_table, chromosome, assembly, "positions written")

	if len(pos_files)==0:
		print("pid %d: no positions found (?)" % os.getpid())
		return

	time1 = time.time()
	print("\t\t\t chromosome %s positions found in %.3f mins "%(chromosome,float(time1-time0)/60))

	# concatenate all input, remove duplicates
	concat_pos = "all_pos.txt"
	subprocess.call(["bash","-c", " cat %s | sort -g | uniq > %s" % (" ".join(pos_files), concat_pos) ])
	time2 = time.time()
	print("\t\t\t chromosome %s positions concatenated to %s in %.3f mins "%(chromosome, concat_pos, float(time2-time1)/60))

	tsv_list = make_location_tsv(cursor, chromosome, concat_pos)
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
	somatic_temp_tables = other_args[1]
	home = os.getcwd()

	for chrom in chromosomes:

		# make a workdir and move there
		workdir  = "chrom_%s" % chrom
		workpath = "{}/locations/{}".format(home,workdir)
		if not os.path.exists(workpath): os.makedirs(workpath)
		os.chdir(workpath)

		time0 = time.time()
		print("====================")
		print("reorganizing locations for chrom %s, pid %d" % (chrom, os.getpid()))
		reorganize_locations(cursor, ref_assembly, chrom, somatic_temp_tables)
		time1 = time.time()
		print(("chrom %s  done in %.3f mins" % (chrom, float(time1-time0)/60)), os.getpid())

	cursor.close()
	db.close()

	return

#########################################
#########################################
#########################################
def main():

	print("Disabled. Loads location tables without checking.")
	exit()

	ref_assembly = 'hg19' # this is the assembly I would like to see for all coords in location tables
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	somatic_temp_tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	chromosomes = [str(i) for i in range(1,13)] + ["Y"] + [str(i) for i in range(22,12,-1)] + ["X"]
	number_of_chunks = 12  # myISAM does not deadlock

	#chromosomes = ["Y"]
	#number_of_chunks = 1
	parallelize(number_of_chunks, reorganize, chromosomes, [ref_assembly, somatic_temp_tables], round_robin=True)


#########################################
if __name__ == '__main__':
	main()
