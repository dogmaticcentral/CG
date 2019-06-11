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

###################################################
def binvals_fit(bin_centerpoints, bin_vals):
	# decimation - Bezier chokes on too many points
	superbin_size = int(len(bin_centerpoints)/20)
	x = [0]
	y = [0]
	for i in range(0, len(bin_centerpoints), superbin_size):
		n = len(bin_centerpoints[i:i+superbin_size])
		x.append(sum(bin_centerpoints[i:i+superbin_size])/n)
		y.append(sum(bin_vals[i:i+superbin_size])/n)

	xnorm = x[-1] # bezier apparently does not like really big numbers (as in 10^5)
	x = [e/xnorm for e in x]

	nodes = np.asfortranarray([x,y])
	curve = bezier.Curve(nodes, degree=2)

	# plt.scatter(x, y)
	# plt.title("test")
	# plt.ylabel('fraction of donors with mutation in sample')
	# plt.xlabel('cdna length')
	# ax = plt.gca()
	# curve.plot(num_pts=256, ax=ax, color = 'red')
	# plt.ylim([0.0, 0.32])
	# p = [float(x[-2]/2), float(161745)/norm]
	# vertical_line = bezier.Curve(np.asfortranarray([[float(p[0]), float(p[0])],[0.0, 0.75]]), degree=1)
	# vertical_line.plot(num_pts=25, ax=ax, color = 'red')
	# vertical_line = bezier.Curve(np.asfortranarray([[float(p[1]), float(p[1])],[0.0, 0.75]]), degree=1)
	# vertical_line.plot(num_pts=25, ax=ax, color = 'blue')
	# plt.show()
	# print(" %d   %.3f" % (p[0]*norm, curve_value_at_x (curve, p[0])) )
	# print(" %d   %.3f" % (p[1]*norm, curve_value_at_x (curve, p[1])) )
	# exit()

	return curve, xnorm

####################################################
def avg_expectations(cursor, stats_id, tumor_short, number_of_donors):
	avg   = [0]
	stdev = [0]
	bin_centerpoints= [0]

	qry  =  "select parameters, stats from stats where stats_id='%s' " % stats_id
	qry += "and parameters like '%s%%' order by id" % tumor_short
	for params, stats in hard_landing_search(cursor,qry):
		cancer, bin_centerpoint = params.split(";")[:2]
		bin_centerpoint = float(bin_centerpoint) + 500
		a, s = stats.split(";")
		bin_centerpoints.append(bin_centerpoint)
		avg.append(float(a)/number_of_donors)
		stdev.append(float(s)/number_of_donors)


	# Bezier curve survives all kinds of numerical idiocy - and does not take initial guess
	# however, it chokes on too many points to fit
	curve_avg, xnorm_avg  = binvals_fit(bin_centerpoints, avg)
	curve_stdev, xnorm_stdev = binvals_fit(bin_centerpoints, stdev)

	return [bin_centerpoints, avg, stdev, curve_avg, xnorm_avg, curve_stdev, xnorm_stdev]


####################################################
def curve_value_at_x (curve, x):
	#  strategy = bezier.curve.IntersectionStrategy.ALGEBRAIC fails with
	#  The only degrees supported at this time are 0, 1, 2, 3 and 4 (degree=11)

	# crashing on arbitrarty intersection points - never figured out why
	# File "src/bezier/_speedup.pyx", line 485, in bezier._speedup.curve_intersections
	# ValueError: Curve intersection failed to converge to approximately linear subdivisions after 20 iterations.
	# solved by renormalizing the x-range to [0, 1.0] interval

	strategy = bezier.curve.IntersectionStrategy.GEOMETRIC
	vertical_line = bezier.Curve(np.asfortranarray([[float(x), float(x)],[-0.75, 0.75]]), degree=1)
	intersections = curve.intersect(vertical_line,strategy=strategy) # why is this plural?  why two parameters?
	s_vals = np.asfortranarray(intersections[0, :])
	point = curve.evaluate_multi(s_vals)
	# point is given as [[x],[y]]

	if point.size==2:
		return point[1][0]
	else:
		return None


####################################################

def reactome_groups_in_tumor(cursor, table, number_of_donors, number_of_genes_mutated,  stats_id, gene_groups, cdna_length, plot=False):

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

	ret = avg_expectations(cursor,  stats_id, tumor_short, number_of_donors)
	[bin_centerpoints, avg, stdev, curve_avg, xnorm_avg, curve_stdev, xnorm_stdev] = ret

	# for all reactome groups, check the coverage in all tables and show the z-score
	red_points_x = []
	red_points_y = []
	green_points_x = []
	green_points_y = []

	for parent, group in gene_groups.items():
		group_size = len(group)
		if group_size==0: continue

		pathway = get_pathway_name(cursor, parent)
		gene_string = ",".join([quotify(g) for g in group])
		qry  = "select count(distinct icgc_sample_id) from %s " % table
		qry += "where pathogenicity_estimate=1 and reliability_estimate=1 "
		qry += "and gene_symbol in (%s)" % gene_string
		group_mutated = error_intolerant_search(cursor, qry)[0][0]
		scaled_donors_affected = float(group_mutated)/number_of_donors

		cdnal = cdna_length[parent]
		avg_interpolated   = curve_value_at_x(curve_avg, cdnal/xnorm_avg)
		if not avg_interpolated: continue # how does that happen?
		stdev_interpolated = curve_value_at_x(curve_stdev, cdnal/xnorm_stdev)
		if not stdev_interpolated: continue
		z = 0
		if stdev_interpolated>0: z = (scaled_donors_affected-avg_interpolated)/stdev_interpolated
		if abs(z)>3.0:
			outf.write("\n{}\n".format(pathway))
			outf.write("\t number of genes:  %3d   combined cdna length: %d \n" % (group_size, cdnal) )
			outf.write("\t number of donors affected:  %3d \n" % group_mutated)
			outf.write("\t  expected donors affected:  %6.2f      stdev:  %6.2f    z:  %6.2f \n" %
					(avg_interpolated*number_of_donors, stdev_interpolated*number_of_donors, z))

			red_points_x.append(float(cdnal)/xnorm_avg)
			red_points_y.append(scaled_donors_affected)
		else:
			green_points_x.append(float(cdnal)/xnorm_avg)
			green_points_y.append(scaled_donors_affected)

	outf.close()

	# inspect using mathplotlib
	if plot:
		# avg
		# plt.ylim(-0.01, 1.0)
		xvals = [c/xnorm_avg for c in bin_centerpoints]
		# the actual groups of genes
		plt.scatter(green_points_x, green_points_y, c='green', alpha=0.8)
		plt.scatter(red_points_x, red_points_y, c='red', alpha=0.8)
		plt.scatter(xvals, avg,  alpha=0.3)
		# gca ("get current axes") helper function:
		ax = plt.gca()
		# usually matplot lib will put the tickmarks at 0, 0,2, 0.4, 0.6, 0.8, 1.0] what if it does not?
		# the trick from stack overflow ist to call
		# fig, ax = plt.subplots()
		# gig.canvas.draw(plt.show(block=False))
		# and then can get labels with ax.get_xticklabels(labels)
		# but htis pops up a window that I do not want
		labels = [int(float(l)*xnorm_avg) for l in [0, 0,2, 0.4, 0.6, 0.8, 1.0]]
		ax.set_xticklabels(labels)
		plt.title("{}, {}, {} - avg".format(tumor_short, number_of_donors, number_of_genes_mutated))
		plt.ylabel('fraction of donors with mutation in sample')
		plt.xlabel("cdna length")
		# need axis object from plt to plot bezier
		curve_avg.plot(num_pts=256, ax=ax, color='red')

		# plt.show()
		plt.savefig("{}/pathways.avg.png".format(outdir))
		plt.clf()

		# stdev
		# plt.ylim(-0.01, 0.15)
		plt.scatter(xvals, stdev)
		ax = plt.gca()
		plt.title("{}, {}, {} - stdev".format(tumor_short, number_of_donors, number_of_genes_mutated))
		plt.ylabel('stdev - fraction of number of patents')
		ax.set_xticklabels(labels)
		plt.xlabel("cdna length")
		#
		curve_stdev.plot(num_pts=256, ax=ax,  color='red')
		plt.savefig("{}/pathways.stdev.png".format(outdir))
		plt.clf() # clear function - otherwise matplotlib keeps piling the dots from prev graphs


	return


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
def plot_binvals(bin_centerpoints, bin_vals, tumor_short, label, vert_line_at_x = None):

	curve, xnorm = binvals_fit(bin_centerpoints, bin_vals)
	xvals = [c/xnorm for c in bin_centerpoints]
	plt.scatter(xvals, bin_vals, marker="x", alpha=0.5)
	plt.title("{}  - {}".format(tumor_short, label))
	plt.ylabel('fraction of donors with mutation in sample')
	plt.xlabel('cdna length')
	ax = plt.gca()
	curve.plot(num_pts=256, ax=ax, color = 'red')

	if vert_line_at_x:
		x = vert_line_at_x/xnorm
		plt.ylim([0.0, 0.3])
		vertical_line = bezier.Curve(np.asfortranarray([[float(x), float(x)],[0.0, 0.75]]), degree=1)
		vertical_line.plot(num_pts=25, ax=ax, color = 'red')

	plt.show()

##########
def plot_sim_results(cursor, stats_id, tumor_short):
	avg   = [0]
	stdev = [0]
	bin_centerpoints= [0]

	qry  =  "select parameters, stats from stats where stats_id='%s' " % stats_id
	qry += "and parameters like '%s%%' order by id" % tumor_short
	for params, stats in hard_landing_search(cursor,qry):
		cancer, bin_centerpoint = params.split(";")[:2]
		bin_centerpoint = float(bin_centerpoint) + 500
		a, s = stats.split(";")
		bin_centerpoints.append(bin_centerpoint)
		avg.append(float(a))
		stdev.append(float(s))

	plot_binvals(bin_centerpoints, avg, tumor_short, "average")
	plot_binvals(bin_centerpoints, stdev, tumor_short, "stdev")


####################################################
def main():

	plot = True

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor, 'icgc')
	stats_id = "RSSCcdna2"

	tables = get_somatic_variant_tables(cursor)
	# tables = ['PRAD_simple_somatic','GACA_simple_somatic','MALY_simple_somatic',
	# 		  'COCA_simple_somatic', 'AML_simple_somatic', 'LMS_simple_somatic',
	# 		  'PACA_simple_somatic', 'BRCA_simple_somatic']
	#tables = ['PRAD_simple_somatic']
	# for table in tables:
	# 	tumor_short = table.split("_")[0]
	# 	plot_sim_results(cursor, stats_id, tumor_short)
	# exit()

	gene_groups = find_gene_groups(cursor)
	# note the plural: groupS
	cdna_length = gene_groups_cdna_length(cursor, gene_groups)

	# for all tables fit the Bezier curve to avg and stdev
	number_of_donors = {}
	number_of_genes_mutated = {}
	for table in tables:
		number_of_donors[table] = len(get_donors(cursor, table))
		number_of_genes_mutated[table] = get_number_of_genes_affected(cursor, table)
		print(table, number_of_donors[table], number_of_genes_mutated [table])

	# the main loop
	for table in tables:
		reactome_groups_in_tumor(cursor, table, number_of_donors[table], number_of_genes_mutated[table],
								stats_id, gene_groups, cdna_length, plot=plot)

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
