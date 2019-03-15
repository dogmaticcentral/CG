#! /usr/bin/python
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
###############################
def protein_mutations (cursor, tables, gene_symbol):

	p53_wt = []
	p53_mut = []
	for table in tables:

		tumor_short = table.split("_")[0]

		qry  = "select g.icgc_mutation_id, m.icgc_specimen_id, m.chromosome "
		qry += "from mutation2gene g,  %s m " % table
		qry += "where g.gene_symbol='%s' " % gene_symbol
		qry += "and g.icgc_mutation_id = m.icgc_mutation_id "
		qry += "and m.pathogenic_estimate=1 and m.reliability_estimate=1"

		ret = search_db(cursor,qry)

		if not ret: continue # no mutations here

		donor_rows = {}
		p53_status_per_specimen = {}
		specimen_seen = {}
		mutations_per_specimen = {}
		for line in ret:
			[mutation_id, specimen_id, chromsome] = line
			qry = "select consequence, aa_mutation from mutations_chrom_%s " % chromsome
			qry += "where icgc_mutation_id='%s' " %  mutation_id
			ret2 = search_db(cursor,qry)
			if not ret2:
				search_db(cursor,qry,verbose=True)
				exit(1)
			[consequence, aa_change] = ret2[0]
			if consequence!='missense': continue

			###################
			# specimen related info
			if not specimen_seen.has_key(specimen_id):
				specimen_seen[specimen_id] = True
				p53_status_per_specimen[specimen_id] = find_53_status(cursor, tumor_short, specimen_id)

			#print specimen_id, "    ", aa_change_cleanup(cursor, aa_change), "    ", p53_status_per_specimen[specimen_id]
			if p53_status_per_specimen[specimen_id][0] == 'pathogenic':
				p53_mut.append(aa_change_cleanup(cursor, aa_change))
			else:
				p53_wt.append(aa_change_cleanup(cursor, aa_change))

	# TODO: note in the paper some mutations recurrent
	p53_wt = set(p53_wt)
	p53_mut = set(p53_mut)
	#return list(p53_wt.difference(p53_mut)), list(p53_mut.difference(p53_wt))
	return list(p53_wt.difference(p53_mut)), list(p53_mut.difference(p53_wt)), list(p53_wt.union(p53_mut))

###############################
def main():



	gene = 'RPL11'
	pdb_file       = "/home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.%s.pdb" % gene
	clustering_prog = "/home/ivana/c-utils/pdb_clust/pc"

	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")

	p53_wt, p53_mut, both = protein_mutations (cursor, tables, gene)


	for background in ['p53_wt', 'p53_mut', 'both']:

		mutations = eval(background)

		# clustering input
		mut_positions_file = "{}.{}.clust.input".format(gene,background)
		outf = open (mut_positions_file,"w")
		for m in mutations:
			outf.write ("%s  %s\n" % (m[0],m[1:-1]))
		outf.close()

		# clustering run
		clustering_output  =  "{}.{}.clust.output".format(gene,background)
		cmd = "{} {} - {} 4.5 > {}".format(clustering_prog, pdb_file, mut_positions_file, clustering_output)
		subprocess.call(["bash","-c", cmd])

		# parse clustering output
		zscore, isolated, clusters = parse_clust_out(clustering_output)
		pymol_file = "{}.{}.clust.pml".format(gene,background)
		basic_pymol_input(pdb_file, isolated, clusters, pymol_file)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
