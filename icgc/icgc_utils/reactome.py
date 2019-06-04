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

# https://networkx.github.io/documentation/stable/index.html
import networkx as nx
from icgc_utils.mysql import  *

def build_reactome_graph(cursor, verbose=False):
	# feed the parent/child pairs as edges into graph
	ret = hard_landing_search(cursor, 'select parent, child from reactome_hierarchy')
	graph = nx.DiGraph(ret) # directed graph
	if verbose:
		print("graph is directed: ", graph.is_directed())
		print("number of edges:", len(graph.edges))
		print("graph is multigraph:", graph.is_multigraph())
	try:
		edges = nx.find_cycle(graph)
	except:
		if verbose: print("hooray, no cycles found")

	return graph


def get_roots(graph):
	# graph.in_degree is a list of pairs, rather than a method
	return [name for name, indegree in graph.in_degree if indegree==0]
