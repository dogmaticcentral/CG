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
def parse_clust_out(clustering_output):
	isolated = []
	clusters = []

	inf = open(clustering_output,"r")
	reading_isolated = False
	reading_cluster = False
	list_of_res = []
	for line in inf:
		line = line.strip()
		if len(line)==0: continue
		if 'z-score' in line:
			zscore = float(line.split()[-1])
			if reading_cluster:
				clusters.append(list_of_res)
		elif 'isolated' in line:
			reading_isolated = True
			reading_cluster = False
			list_of_res = []
			continue
		elif  'cluster size' in line:
			if reading_isolated:
				isolated = list_of_res
			elif reading_cluster:
				clusters.append(list_of_res)
			reading_isolated = False
			reading_cluster = True
			list_of_res = []
			continue
		elif reading_cluster or reading_isolated:
			list_of_res.append(line)

	clusters = sorted(clusters, key=lambda cluster: -len(cluster))

	return zscore, isolated, clusters

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
		subprocess.call(["bash","-c", cmd])

		# parse clustering output
		zscore, isolated, clusters = parse_clust_out(clustering_output)

		# pymol input
		pymol_file = "{}.{}.clust.pml".format(gene,background)
		outf = open(pymol_file,"w")

		outf.write("\n")
		outf.write("load %s\n" %pdbf_file)
		outf.write("\nbg_color  white\nhide all\nshow cartoon\ncolor white\n")


		outf.write("\n")
		outf.write("#isolated\n")
		outf.write("select isolated, resi %s\n" % ("+".join(isolated)) )
		outf.write("show spheres, isolated\n")

		outf.write("color gray, isolated\n")

		outf.write("\n")
		outf.write("#clusters\n")
		outf.write("\n")
		outf.write("select red_clust, resi %s\n" % ("+".join(clusters[0])) )
		outf.write("color red, red_clust\n")
		outf.write("show spheres, red_clust\n")
		outf.write("\n")
		outf.write("select blue_clust, resi %s\n" % ("+".join(clusters[1])) )
		outf.write("color blue, blue_clust\n")
		outf.write("show spheres, blue_clust\n")
		outf.write("\n")
		outf.write("select green_clust, resi %s\n" % ("+".join(clusters[2])) )
		outf.write("color green, green_clust\n")
		outf.write("show spheres, green_clust\n")

		for cluster in clusters[3:]:
			outf.write("\n")
			clust_id = "cluster_%d" % clusters.index(cluster)
			outf.write("select %s, resi %s\n" % (clust_id, "+".join(cluster)) )
			outf.write("color orange, %s\n" % clust_id)
			outf.write("show spheres, %s\n" % clust_id)

			outf.write("\ndeselect\n")


		outf.close()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
