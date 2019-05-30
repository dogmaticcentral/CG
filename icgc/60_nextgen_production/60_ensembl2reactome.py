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

# awk -F '\t' '$1~"ENSG" &&  4!="IEA"  && $6=="Homo sapiens" \
# {printf "%s\t%s\t%s\n", $1, $2, $4}' Ensembl2Reactome.txt  > Ensembl2Reactome.human.txt

# some 'pathways' are not really pathways ("S45 mutants of beta-catenin aren't phosphorylated",
# "yoclonic epilepsy of Lafora", )
# many are description of disregulated events ("Defective ... causes ...")
# several describe HIV related events (""), adn other disease scenarios
# ("Influenza Virus Induced Apoptosis", )
# there is a Complete list of pathways and Pathways hierarchy relationship
# build tree, select pathways we are interested in

import time

from icgc_utils.common_queries import *
from icgc_utils.processes import *
from config import Config
from random import sample
from math import sqrt

def make_tables(cursor):
	# def reactome pathway names varchar(20) vs text
	# def ensembl 2 pathway varchar(20) vs varchar(20)
	return

####################################################
def main():


	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
