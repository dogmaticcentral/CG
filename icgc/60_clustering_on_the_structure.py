#! /usr/bin/python
import subprocess


from icgc_utils.common_queries import *


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
			if p53_status_per_specimen[specimen_id][0] == 'wt':
				p53_wt.append(aa_change_cleanup(cursor, aa_change))
			else:
				p53_mut.append(aa_change_cleanup(cursor, aa_change))

	# TODO: note in the paper some mutations recurrent
	p53_wt = set(p53_wt)
	p53_mut = set(p53_mut)
	#return list(p53_wt.difference(p53_mut)), list(p53_mut.difference(p53_wt))
	return list(p53_wt.difference(p53_mut)), list(p53_mut.difference(p53_wt)), list(p53_wt.union(p53_mut))




###############################
def main():

	gene = 'RPL5'
	pdbf_file       = "/home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.%s.pdb" % gene
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
		cmd = "{} {} - {} 4.5 > {}".format(clustering_prog, pdbf_file, mut_positions_file, clustering_output)
		print cmd
		subprocess.call(["bash","-c", cmd])

		# parse clustering output

		# pymol input



	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
