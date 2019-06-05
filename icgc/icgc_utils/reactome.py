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
from icgc_utils.common_queries import *

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


def count_successors(graph, pthwy_id):
	return len(list(graph.successors(pthwy_id)))


################
def genes_in_subgraph(cursor, graph, parent_id):
	genes = []
	# this is the whole subtree
	descendants = [pid for pid in nx.dfs_preorder_nodes(graph, parent_id) if count_successors(graph, pid) == 0]

	desc_id_string = ",".join([quotify(d) for d in descendants])
	qry  = "select distinct(gene_symbol) from  hgnc2reactome "
	qry += "where reactome_pathway_id in (%s)" % desc_id_string
	ret = error_intolerant_search(cursor,qry)
	if not ret:
		print("possible problem in Reactome: no associated genes found for ", desc_id_string)
		return []
	return [r[0] for r in ret]


################
def ensembl_genes_in_subgraph(cursor, graph, parent_id):
	genes = []
	# this is the whole subtree
	descendants = [pid for pid in nx.dfs_preorder_nodes(graph, parent_id) if count_successors(graph, pid) == 0]

	desc_id_string = ",".join([quotify(d) for d in descendants])
	qry = "select distinct(ensembl_gene_id) from  ensembl2reactome "
	qry += "where reactome_pathway_id in (%s)" % desc_id_string
	ret = error_intolerant_search(cursor,qry)
	if not ret:
		print("possible problem in Reactome: no associated genes found for ", desc_id_string)
		return []
	return [r[0] for r in ret]


def get_pathway_name(cursor, node_id):
	qry  = "select name from reactome_pathways "
	qry += "where reactome_pathway_id = '%s'" % node_id
	ret = error_intolerant_search(cursor,qry)
	if not ret: return "name not found"
	return ret[0][0]

def get_pathway_names(cursor, list_of_node_ids):
	qry  = "select * from reactome_pathways "
	qry += "where reactome_pathway_id in (%s)" % ",".join([quotify(i) for i in list_of_node_ids])
	ret = error_intolerant_search(cursor,qry)
	if not ret: return {}
	return dict(ret)

