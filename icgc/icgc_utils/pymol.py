
def basic_pymol_input(pdb_file, isolated, clusters, output_file):
	# pymol input
	outf = open(output_file,"w")

	outf.write("\n")
	outf.write("load %s\n" %pdb_file)
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
