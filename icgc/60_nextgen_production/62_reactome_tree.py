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
from icgc_utils.common_queries import quotify
from icgc_utils.mysql import  *
from config import Config
# https://networkx.github.io/documentation/stable/index.html
import networkx as nx

####################################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor, 'icgc')

	# are there children with multiple parents? Yes. So I need some kind of
	# directed graph, rather tha a tree.
	qry = "select child, count(distinct parent) as ct from reactome_hierarchy "
	qry += "group by child having ct>1"
	ret = search_db(cursor, qry)
	print("number of children with multiple parents:", len(ret))

	# feed the parent/child pairs as edges into graph
	ret = hard_landing_search(cursor, 'select parent, child from reactome_hierarchy')
	graph = nx.DiGraph(ret) # directed graph
	print("graph is directed: ", graph.is_directed())
	print("number of edges:", len(graph.edges))
	print("graph is multigraph:", graph.is_multigraph())
	try:
		edges = nx.find_cycle(graph)
	except:
		print("hooray, no cycles found")

	# graph.in_degree is a list of pairs, rather than a method

	# candidate roots
	zero_in_degee_nodes = [name for name, indegree in graph.in_degree if indegree==0]

	node_id_string = ",".join([quotify(z) for z in zero_in_degee_nodes])
	qry_template = "select * from reactome_pathways where reactome_pathway_id in (%s)"
	root_names =  hard_landing_search(cursor, qry_template% node_id_string)
	print("zero in-degree nodes:")
	for pthwy_id, name in root_names:
		print(pthwy_id, name)
		# this is the whole subtree
		# children = [node for node in nx.dfs_preorder_nodes(graph, pthwy_id)]
		# A successor of n is a node m such that there exists a directed edge from n to m.
		children = [node for node in graph.successors(pthwy_id)]
		node_id_string = ",".join([quotify(z) for z in children])
		children_names = hard_landing_search(cursor, qry_template% node_id_string)
		for child_id, child_name in children_names:
			print("\t", child_id, child_name)


	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
