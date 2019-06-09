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

from config import Config
from icgc_utils.reactome import *
from icgc_utils.common_queries import *
import numpy as np
from matplotlib import pyplot as plt
import bezier

####################################################
def atomic_groups(cursor, graph, root, groups):
	children = [node for node in graph.successors(root)]
	if len(children)==0: return False
	for child in children:
		# all genes that this sub-graph encompasses
		genes = list(filter(lambda g: g!='TP53', genes_in_subgraph(cursor, graph, child)))
		# if they are fewer than 100, we call that an atomic group
		if len(genes)<100:
			groups[child] = genes
			continue
		# otherwise we keep subdividing
		if not atomic_groups(cursor, graph, child, groups): # no further subdivisions
			groups[child] = genes
			continue

	return groups

####################################################
def avg_expectations(cursor, tumor_short, number_of_donors):
	avg   = [0]
	stdev = [0]
	bin_centerpoints= [0]

	qry = "select parameters, stats from stats where stats_id='RSSCcdna' and parameters like '%s%%'" % tumor_short
	for params, stats in hard_landing_search(cursor,qry):
		cancer, bin_centerpoint = params.split(";")
		bin_centerpoint = float(bin_centerpoint)+ 500
		a, s = stats.split(";")
		bin_centerpoints.append(bin_centerpoint)
		avg.append(float(a)/number_of_donors)
		stdev.append(float(s)/number_of_donors)

	# Bezier curve survives all kinds of numerical idiocy - and does not take initial guess
	# however, it chokes on too many points to fit
	nodes = np.asfortranarray([bin_centerpoints, avg])
	curve_avg = bezier.Curve(nodes, degree=2)
	nodes = np.asfortranarray([bin_centerpoints, stdev])
	curve_stdev = bezier.Curve(nodes, degree=2)

	return [bin_centerpoints, avg, stdev, curve_avg, curve_stdev]


####################################################
def curve_value_at_x (curve, x):
	vertical_line = bezier.Curve(np.asfortranarray([[float(x), float(x)],[0.0, 0.75]]), degree=1)
	intersections = curve.intersect(vertical_line) # why is this plural?  why two parameters?
	s_vals = np.asfortranarray(intersections[0, :])
	point = curve.evaluate_multi(s_vals)
	# point is given as [[x],[y]]

	if point.size==2:
		return point[1][0]
	else:
		return None


####################################################
def reactome_groups_in_tumor(cursor, table, number_of_donors, number_of_genes_mutated,  gene_groups, cdna_length, plot=False):

	if not number_of_donors:
		print("no samples for %s (?)"% table)
		return
	tumor_short = table.split("_")[0]
	outdir = "gene_groups/{}".format(tumor_short)
	if not os.path.exists(outdir): os.makedirs(outdir)
	outf = open("{}/pathways.txt".format(outdir), "w")

	print("\n\n=====================")
	print("{}     number of donors: {}    number of genes mutated: {}  ".
			format(tumor_short, number_of_donors, number_of_genes_mutated))
	outf.write("{}     number of donors: {}    number of genes mutated: {}  \n".
			format(tumor_short, number_of_donors, number_of_genes_mutated))

	[bin_centerpoints, avg, stdev, curve_avg, curve_stdev] = avg_expectations(cursor, tumor_short, number_of_donors)

	# for all reactome groups, check the coverage in all tables and show the z-score
	red_points_x = []
	red_points_y = []
	green_points_x = []
	green_points_y = []

	for parent, group in gene_groups.items():
		group_size = len(group)
		if group_size==0: continue
		if cdna_length[parent]>=bin_centerpoints[-1]: continue

		pathway = get_pathway_name(cursor, parent)
		gene_string = ",".join([quotify(g) for g in group])
		qry  = "select count(distinct icgc_sample_id) from %s " % table
		qry += "where pathogenicity_estimate=1 and reliability_estimate=1 "
		qry += "and gene_symbol in (%s)" % gene_string
		group_mutated = error_intolerant_search(cursor, qry)[0][0]
		scaled_donors_affected = float(group_mutated)/number_of_donors

		avg_interpolated   = curve_value_at_x(curve_avg, cdna_length[parent])
		if not avg_interpolated: continue # how does that happen?
		stdev_interpolated = curve_value_at_x(curve_stdev, cdna_length[parent])
		if not stdev_interpolated: continue
		z = 0
		if stdev_interpolated>0: z = (scaled_donors_affected-avg_interpolated)/stdev_interpolated
		if abs(z)>3.0:
			outf.write("\n{}\n".format(pathway))
			outf.write("\t number of genes:  %3d   combined cdna length: %d \n" % (group_size, cdna_length[parent]) )
			outf.write("\t number of donors affected:  %3d \n" % group_mutated)
			outf.write("\t  expected donors affected:  %6.2f      stdev:  %6.2f    z:  %6.2f \n" %
					(avg_interpolated*number_of_donors, stdev_interpolated*number_of_donors, z))

			red_points_x.append(group_size)
			red_points_y.append(scaled_donors_affected)
		else:
			green_points_x.append(group_size)
			green_points_y.append(scaled_donors_affected)

	outf.close()

	# inspect using mathplotlib
	if plot:
		# avg
		# plt.ylim(-0.01, 1.0)
		plt.scatter(bin_centerpoints, avg)
		plt.title("{}, {}, {} - avg".format(tumor_short, number_of_donors, number_of_genes_mutated))
		plt.ylabel('fraction of donors with mutation in sample')
		plt.xlabel('sample size')
		# need axis object from plt to plot bezier
		# gca ("get current axes") helper function:
		ax = plt.gca()
		curve_avg.plot(num_pts=256, ax=ax)

		# the actual groups of genes
		plt.scatter(red_points_x, red_points_y, c='red')
		plt.scatter(green_points_x, green_points_y, c='green')
		plt.savefig("{}/pathways.png".format(outdir))
		plt.clf()
		# stdev
		# plt.ylim(-0.01, 0.15)
		# plt.scatter(sample_sizes, stdev)
		# plt.title("{}, {}, {} - stdev".format(tumor_short, number_of_donors, number_of_genes_mutated))
		# plt.ylabel('fraction of donors with mutation in sample')
		# # need axis objet from plt to plot bezier
		# # gca ("get current axes") helper function:
		# ax = plt.gca()
		# curve_stdev.plot(num_pts=256, ax=ax)
		# plt.show()

	return


####################################################
def find_gene_groups(cursor):
	# feed the parent/child pairs as edges into graph
	graph = build_reactome_graph(cursor, verbose=True)
	# candidate roots
	zero_in_degee_nodes = get_roots(graph)
	root_name = get_pathway_names(cursor, zero_in_degee_nodes)

	gene_groups = {}
	for root in zero_in_degee_nodes:
		if 'disease' in root_name[root].lower(): continue
		atomic_groups(cursor, graph, root, gene_groups)

	return gene_groups

####################################################
def gene_groups_cdna_length(cursor, gene_groups):
	cdna_length = {}
	for parent, group in gene_groups.items():

		if len(group)==0:
			cdna_length[parent]=0
			continue
		# note the singular - defined in common queries
		total_length = gene_group_cdna_length(cursor, group)

		cdna_length[parent] = total_length

	return cdna_length


####################################################
from lmfit import Minimizer, Parameters, report_fit


def plot_sim_results(cursor, tumor_short):
	avg   = [0]
	stdev = [0]
	bin_centerpoints= [0]

	qry = "select parameters, stats from stats where stats_id='RSSCcdna' and parameters like '%s%%' order by id" % tumor_short
	for params, stats in hard_landing_search(cursor,qry):
		cancer, bin_centerpoint = params.split(";")
		bin_centerpoint = float(bin_centerpoint) + 500
		a, s = stats.split(";")
		bin_centerpoints.append(bin_centerpoint)
		avg.append(float(a))
		stdev.append(float(s))

	for i in range(30):
		print(i, bin_centerpoints[i], avg[i])

	# decimation - Bezier chokes on too many points
	x = [0]
	y = [0]
	for i in range(0, len(bin_centerpoints), 10):
		n = len(bin_centerpoints[i:i+10])
		x.append(sum(bin_centerpoints[i:i+10])/n)
		y.append(sum(avg[i:i+10])/n)

	nodes = np.asfortranarray([x,y])
	curve_avg = bezier.Curve(nodes, degree=2)

	plt.scatter(bin_centerpoints, avg)
	plt.title("{}  - avg".format(tumor_short))
	plt.ylabel('fraction of donors with mutation in sample')
	plt.xlabel('cdna length')
	ax = plt.gca()
	curve_avg.plot(num_pts=256, ax=ax, color = 'red')

	plt.show()

####################################################
def main():

	plot = True

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor, 'icgc')


	tables = get_somatic_variant_tables(cursor)
	for table in tables:
		tumor_short = table.split("_")[0]
		if tumor_short in ['MELA', 'ESAD']: continue
		plot_sim_results(cursor, tumor_short)
	exit()


	gene_groups = find_gene_groups(cursor)
	# note the plural: groupS
	cdna_length = gene_groups_cdna_length(cursor, gene_groups)

	# for all tables fit the Bezier curve to avg and stdev
	number_of_donors = {}
	number_of_genes_mutated = {}
	for table in tables:
		number_of_donors[table] = len(get_donors(cursor, table))
		number_of_genes_mutated [table] = get_number_of_genes_affected(cursor, table)
		print(table, number_of_donors[table], number_of_genes_mutated [table])

	# the main loop
	for table in tables:
		reactome_groups_in_tumor(cursor, table,
								 number_of_donors[table], number_of_genes_mutated[table],
								 gene_groups, cdna_length, plot=plot)

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
