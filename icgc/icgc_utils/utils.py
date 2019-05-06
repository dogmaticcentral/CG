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


#########################################
def find_clusters(string_pairs):
	clusters = []
	for pair in string_pairs:
		new_set = set(pair)
		placed = False
		for cluster in clusters:
			# if intersection => join and break
			if len(cluster.intersection(new_set))==0: continue
			cluster |= new_set
			placed = True
			break
		if not placed:
			clusters.append(new_set)
	return clusters
