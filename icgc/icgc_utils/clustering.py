###############################
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

