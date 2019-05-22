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

from icgc_utils.common_queries import *
from icgc_utils.pymol import *
from icgc_utils.clustering import *
from config import Config

###############################
def main():


	if len(sys.argv)<3:
		print("usage: %s <gene symbol>  <pdb file> " % sys.argv[0])
		exit()

	gene    = sys.argv[1].upper()
	pdb_file = sys.argv[2]

	clustering_prog = "/home/ivana/c-utils/pdb_clust/pc"

	for fnm in [pdb_file, clustering_prog]:
		if os.path.exists(fnm) and os.path.getsize(fnm)>0: continue
		print(fnm,"not found or empty")
		exit(1)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")

	mutations = gnomad_mutations (cursor, gene)

	# clustering input
	mut_positions_file = "{}.{}.clust.input".format(gene,"gnomad")
	outf = open (mut_positions_file,"w")
	for m in mutations:
		outf.write ("%s  %s\n" % (m[0],m[1:-1]))
	outf.close()

	# clustering run
	clustering_output  =  "{}.{}.clust.output".format(gene,"gnomad")
	cmd = "{} {} - {} 4.5 > {}".format(clustering_prog, pdb_file, mut_positions_file, clustering_output)
	print(cmd)
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
