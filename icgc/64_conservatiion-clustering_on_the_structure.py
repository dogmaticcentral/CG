#! /usr/bin/python
import subprocess


from icgc_utils.common_queries import *
from icgc_utils.pymol import *
from icgc_utils.clustering import *

###############################
def main():

	gene = 'RPL11'
	pdb_file        = "/home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.%s.pdb" % gene
	clustering_prog = "/home/ivana/c-utils/pdb_clust/pc"
	conservation_score_file  = "/home/ivana/Dropbox/Sinisa/ribosomal/data/%s/" % gene.lower()
	conservation_score_file += "conservation_animals_plants_fungi/specs_out.score"

	for fnm in [pdb_file,clustering_prog,conservation_score_file]:
		if not os.path.exists(fnm):
			print fnm, "not found"
			exit()

	# clustering input
	mut_positions_file = "{}.{}.clust.input".format(gene,"conservation")
	inf  = open(conservation_score_file,"r")
	outf = open (mut_positions_file,"w")
	for line in inf:
		if line[0]=='%': continue
		fields = line.strip().split()
		if float(fields[2]) <= 0.05:
			outf.write ("%s  %s\n" % (fields[3], fields[4]))
	outf.close()

	# clustering run
	clustering_output  =  "{}.{}.clust.output".format(gene,"conservation")
	cmd = "{} {} - {} 4.5 > {}".format(clustering_prog, pdb_file, mut_positions_file, clustering_output)
	subprocess.call(["bash","-c", cmd])

	# parse clustering output
	zscore, isolated, clusters = parse_clust_out(clustering_output)
	pymol_file = "{}.{}.clust.pml".format(gene,"conservation")
	basic_pymol_input(pdb_file, isolated, clusters, pymol_file)



#########################################
if __name__ == '__main__':
	main()
