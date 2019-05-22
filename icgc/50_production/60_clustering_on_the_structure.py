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

from icgc_utils.common_queries import *
from icgc_utils.pymol import *
from icgc_utils.clustering import *
from config import Config

###############################
def protein_mutations (cursor, tables, gene_symbol, bg_gene):

	canonical_transcript  = {}
	for g in [gene_symbol, bg_gene]:
		canonical_transcript[g] = approved_symbol2ensembl_canonical_transcript(cursor,g)
		if not canonical_transcript[g]:
			print ("canonical transcript not found for", g)
			exit()

	bg_wt = []
	bg_mut = []
	for table in tables:

		tumor_short = table.split("_")[0]

		qry  = "select g.icgc_mutation_id, m.icgc_specimen_id, m.chromosome "
		qry += "from mutation2gene g,  %s m " % table
		qry += "where g.gene_symbol='%s' " % gene_symbol
		qry += "and g.icgc_mutation_id = m.icgc_mutation_id "
		qry += "and m.pathogenicity_estimate=1 and m.reliability_estimate=1"

		ret = error_intolerant_search(cursor,qry)

		if not ret: continue # no mutations here

		bg_status_per_specimen = {}
		specimen_seen = {}
		bgct = canonical_transcript[bg_gene]
		for line in ret:
			[mutation_id, specimen_id, chromosome] = line
			consequence, aa_change  = get_consequence(cursor, chromosome, mutation_id)
			if not 'missense' in consequence: continue
			aa_change = consequence_cleanup(canonical_transcript[gene_symbol], aa_change)
			if ":" in aa_change: continue # mutation in an alternative transcript only
			###################
			# specimen related info
			if not specimen_id  in specimen_seen:
				specimen_seen[specimen_id] = True
				bg_status_per_specimen[specimen_id] = find_background_status(cursor, tumor_short, specimen_id, bg_gene, bgct)
			#print(tumor_short , specimen_id, "  ",  aa_change, "  ", bg_status_per_specimen[specimen_id])
			if bg_status_per_specimen[specimen_id][0] == 'pathogenic':
				bg_mut.append(aa_change)
			else:
				bg_wt.append(aa_change)

	bg_wt = set(bg_wt)
	bg_mut = set(bg_mut)

	return list(bg_wt.difference(bg_mut)), list(bg_mut.difference(bg_wt)), list(bg_wt.union(bg_mut))

###############################
def main():

	if len(sys.argv)<3:
		print("usage: %s <gene symbol>  <pdb file> [<background gene symbol>]" % sys.argv[0])
		exit()
		
	gene    = sys.argv[1].upper()
	pdb_file = sys.argv[2]
	bg_gene = sys.argv[3] if len(sys.argv)>2 else 'TP53'
	
	clustering_prog = "/home/ivana/c-utils/pdb_clust/pc"
	
	for fnm in [pdb_file, clustering_prog]:
		if os.path.exists(fnm) and os.path.getsize(fnm)>0: continue
		print(fnm,"not found or empty")
		exit(1)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#########################
	# which simple somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	#########################
	switch_to_db(cursor,"icgc")

	bg_wt, bg_mut, both = protein_mutations (cursor, tables, gene, bg_gene)

	for background in ['bg_wt', 'bg_mut', 'both']:

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
