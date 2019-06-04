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
from icgc_utils.mysql import *
from icgc_utils.common_queries import *
import numpy as np
from matplotlib import pyplot as plt
import bezier

# todo: bezier for random selections:
# https://pypi.org/project/bezier/, https://github.com/dhermes/bezier/blob/master/src/bezier/curve.py
####################################################
def main():

	plot = True

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor, 'icgc')

	# for all tables fit the Bezier curve to avg and stdev
	tables = get_somatic_variant_tables(cursor)
	number_of_donors = {}
	number_of_genes_mutated = {}
	for table in tables:
		number_of_donors[table] = len(get_donors(cursor, table))
		number_of_genes_mutated [table] = get_number_of_genes_affected(cursor, table)
		print(table, number_of_donors[table], number_of_genes_mutated [table])
		# CMDI avgs behave in a bizarre way, investigate at some other point
		if table=="CMDI": continue

	curve_avg = {}
	curve_stdev = {}
	for table in tables:
		# CMDI avgs behave in a bizarre way, investigate at some other point
		if not number_of_donors[table] :
			print("no samples for %s (?)"% table)
			continue
		tumor_short = table.split("_")[0]
		avg   = []
		stdev = []
		sample_sizes= []
		qry = "select parameters, stats from stats where stats_id='RSSC' and parameters like '%s%%'" % tumor_short
		for params, stats in hard_landing_search(cursor,qry):
			cancer, sample_size = params.split(";")
			a, s = stats.split(";")
			sample_sizes.append(float(sample_size))
			avg.append(float(a)/number_of_donors[table])
			stdev.append(float(s)/number_of_donors[table])

		# fit Bezier curve - it survives all kinds of numerical idiocy
		nodes = np.asfortranarray([sample_sizes, avg])
		curve_avg[table] = bezier.Curve(nodes, degree=2)
		nodes = np.asfortranarray([sample_sizes, stdev])
		curve_stdev[table] = bezier.Curve(nodes, degree=2)

		# inspect using mathplotlib
		if plot:
			# avg
			plt.ylim(-0.01, 1.0)
			plt.scatter(sample_sizes, avg)
			plt.title("{}, {}, {} - avg".format(tumor_short, number_of_donors[table], number_of_genes_mutated[table]))
			plt.ylabel('fraction of donors with mutation in sample')
			# need axis objet from plt to plot bezier
			# gca ("get current axes") helper function:
			ax = plt.gca()
			curve_avg[table].plot(num_pts=256, ax=ax)
			plt.show()

			# stdev
			plt.ylim(-0.01, 0.15)
			plt.scatter(sample_sizes, stdev)
			plt.title("{}, {}, {} - stdev".format(tumor_short, number_of_donors[table], number_of_genes_mutated[table]))
			plt.ylabel('fraction of donors with mutation in sample')
			# need axis objet from plt to plot bezier
			# gca ("get current axes") helper function:
			ax = plt.gca()
			curve_stdev[table].plot(num_pts=256, ax=ax)
			plt.show()


	# for all reactome groups, check the coverage in all tables and show the z-score

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
