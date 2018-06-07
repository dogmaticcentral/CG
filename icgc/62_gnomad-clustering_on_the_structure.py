#! /usr/bin/python
import subprocess


from icgc_utils.common_queries import *
from icgc_utils.pymol import *
from icgc_utils.clustering import *

###############################
def main():

	gene = 'RPL5'
	pdb_file       = "/home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.%s.pdb" % gene
	clustering_prog = "/home/ivana/c-utils/pdb_clust/pc"

	db     = connect_to_mysql()
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")

	mutations = gnomad_mutations (cursor, gene)

	# clustering input
	mut_positions_file = "{}.{}.clust.input".format(gene,"gnomad")
	outf = open (mut_positions_file,"w")
	for m in mutations:
		if m =='M200' : m='L200'
		outf.write ("%s  %s\n" % (m[0],m[1:-1]))
	outf.close()

	# clustering run
	clustering_output  =  "{}.{}.clust.output".format(gene,"gnomad")
	cmd = "{} {} - {} 4.5 > {}".format(clustering_prog, pdb_file, mut_positions_file, clustering_output)
	subprocess.call(["bash","-c", cmd])

	# parse clustering output
	zscore, isolated, clusters = parse_clust_out(clustering_output)
	pymol_file = "{}.{}.clust.pml".format(gene,"gnomad")
	basic_pymol_input(pdb_file, isolated, clusters, pymol_file)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
