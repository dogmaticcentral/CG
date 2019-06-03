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

# https://reactome.org/download-data
# Lowest level pathway diagram / Subset of the pathway
# ENSEMBL to pathways --> save as Ensembl2Reactome.txt
# fields
#     Source database identifier, e.g. UniProt, ENSEMBL, NCBI Gene or ChEBI identifier
#     Reactome Pathway Stable identifier
#     URL
#     Event (Pathway or Reaction) Name
#     Evidence Code
#     Species
# evidence code: for human only TAS (traceable author statement) and
# IEA (inferred electronic annotation) - for most pathways with evidence code IEA, TAS also exists


# some 'pathways' are not really pathways ("S45 mutants of beta-catenin aren't phosphorylated",
# "myoclonic epilepsy of Lafora", )
# many are description of disregulated events ("Defective ... causes ...")
# several describe HIV related events (""), adn other disease scenarios
# ("Influenza Virus Induced Apoptosis", )
# there is a Complete list of pathways and Pathways hierarchy relationship
# build tree, select pathways we are interested in


# awk -F '\t' '$1~"ENSG" &&  4!="IEA"  && $6=="Homo sapiens" \
#  {ct++; printf "%d\t%s\t%s\n", ct, $1, $2}' Ensembl2Reactome.txt  > Ensembl2Reactome.human.tsv

# awk -F '\t' '$1~"R-HSA" {printf "%s\t%s\n", $1, $2}' ReactomePathways.txt > ReactomePathways.human.tsv

# awk -F '\t' '$1~"R-HSA" { ct++; printf "%d\t%s\t%s\n", ct, $1, $2}' \
#   ReactomePathwaysRelation.txt > ReactomePathwaysRelation.human.tsv

# bash> mv Ensembl2Reactome.human.tsv ensembl2reactome.tsv
# bash> mv ReactomePathways.human.tsv reactome_pathways.tsv
# bash> mv ReactomePathwaysRelation.human.tsv reactome_hierarchy.tsv
# bash> mkdir tsvs && mv *.tsv tsvs

import time

from icgc_utils.common_queries import *
from icgc_utils.processes import *
from config import Config
from random import sample
from math import sqrt

def make_tables(cursor, db_name):
	switch_to_db(cursor, db_name)
	# def reactome pathway names varchar(20) vs text
	table =  'reactome_pathways'
	if not check_table_exists(cursor, db_name, table):
		qry = ""
		qry += "  CREATE TABLE  %s (" % table
		qry += "  	 reactome_pathway_id VARCHAR (20) NOT NULL, "
		qry += "     name TEXT NOT NULL, "
		qry += "	 PRIMARY KEY (reactome_pathway_id) "
		qry += ") ENGINE=MyISAM"
		search_db(cursor, qry, verbose=True)

	# def ensembl 2 reactome_pathway varchar(20) vs varchar(20)
	table =  'ensembl2reactome'
	if not check_table_exists(cursor, db_name, table):
		qry = ""
		qry += "  CREATE TABLE  %s (" % table
		qry += "     id INT NOT NULL, "
		qry += "  	 ensembl_gene_id VARCHAR (20) NOT NULL, "
		qry += "  	 reactome_pathway_id VARCHAR (20) NOT NULL, "
		qry += "	 PRIMARY KEY (id) "
		qry += ") ENGINE=MyISAM"
		search_db(cursor, qry, verbose=True)

	# reactome hierarchy
	table = 'reactome_hierarchy'
	if not check_table_exists(cursor, db_name, table):
		qry = ""
		qry += "  CREATE TABLE  %s (" % table
		qry += "     id INT NOT NULL, "
		qry += "  	 parent VARCHAR (20) NOT NULL, "
		qry += "  	 child VARCHAR (20) NOT NULL, "
		qry += "	 PRIMARY KEY (id) "
		qry += ") ENGINE=MyISAM"
		search_db(cursor, qry, verbose=True)

	return

####################################################
def main():


	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	# does nothing if tables exist
	make_tables(cursor, 'icgc')

	for table in ['reactome_pathways', 'ensembl2reactome','reactome_hierarchy']:
		qry = "load data local infile 'tsvs/%s.tsv' into table %s" % (table,table)
		search_db(cursor,qry,verbose=True)


	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
