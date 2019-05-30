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
	print("graph is multigraph: ", graph.is_multigraph())
	try:
		edges = nx.find_cycle(graph)
	except:
		print("hooray, no cycles found")
	

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
