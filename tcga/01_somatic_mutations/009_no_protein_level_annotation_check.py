#!/usr/bin/python

#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2018 Ivana Mihalek.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#



import MySQLdb
from sets import Set
from   tcga_utils.mysql   import  *
import commands


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names = ["ACC", "BLCA", "BRCA", "CESC", "CHOL",  "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",  "MESO", "OV",   "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    sample_type = "metastatic"

    if sample_type == "primary":
        table_name = 'somatic_mutations'
    elif sample_type == "metastatic":
        table_name = 'metastatic_mutations'
    else:
        print "I don't know how to hadndle ", sample_type, " sample types"
        exit(1) # unknown sample type

    # how many cases do we have affected?
    total_cases = 0
    for db_name in db_names:
        print " ** ", db_name
        switch_to_db (cursor, db_name)
        per_db_cases = 0
        qry = "select count(1) from %s where variant_classification='missense_mutation' " % table_name
        qry += " and (aa_change is null or aa_change='')";
        rows = search_db (cursor, qry)
        if rows and rows[0][0] != 0:
            print ">>>>>> ",  rows[0][0]
            per_db_cases = rows[0][0]

        out_of = 0
        qry = "select count(1) from %s where variant_classification='missense_mutation' " % table_name
        rows = search_db (cursor, qry)
        if rows and rows[0][0] != 0:
            out_of = rows[0][0]

        print "number of missing protein annotation cases: ", per_db_cases, "out of", out_of
        print
    
        total_cases += per_db_cases


    print
    print "total number of missing annotation cases: ", total_cases
    print

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
