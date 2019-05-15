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

import os
from icgc_utils.CrossMap import *

#########################################
def find_clusters(string_pairs):
	clusters = []
	for pair in string_pairs:
		new_set = set(pair)
		placed = False
		for cluster in clusters:
			# if intersection => join and break
			if len(cluster.intersection(new_set))==0: continue
			cluster |= new_set
			placed = True
			break
		if not placed:
			clusters.append(new_set)
	return clusters


#########################################
def translate_positions(positions, chromosome, from_assembly, to_assembly, rootname=None):

	if from_assembly == to_assembly:
		return positions
	# GRCh37 and hg19 only differ for MT
	if (from_assembly.lower() in ['grch37', 'hg19']) and (to_assembly.lower() in ['grch37', 'hg19']) and (chromosome != "MT"):
		return positions
	if "grch" in from_assembly.lower():  from_assembly = from_assembly.lower().replace("grch", "GRCh")
	if "grch" in to_assembly.lower():  to_assembly = to_assembly.lower().replace("grch", "GRCh")

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
		new_positions = [int(line.split("\t")[1]) for line in inf.read().split("\n") if len(line.replace(" ",""))>0]

	if len(new_positions) != len(positions):
		print("translation mismatch, {}, chromosome {}, {} to {}".format(os.getcwd(), chromosome, from_assembly, to_assembly))
		#print("process", os.getpid(), "exiting")
		return None
	else:
		# remove aux files
		os.remove(outfile)
		os.remove(outfile_translated)
		os.remove(outfile_translated+".unmap") # this file should probably checked - it should be empty

	return new_positions

#########################################