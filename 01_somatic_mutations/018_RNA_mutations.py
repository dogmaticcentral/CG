#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2016 Ivana Mihalek.
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
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *



#########################################
def main():


    full_name = read_cancer_names ()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG",  "LIHC", "LUAD", "LUSC", "OV",  "PAAD", "PCPG", "PRAD", "REA",
                 "SARC","SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    db_names = ["LUAD"]

   
    table = 'somatic_mutations'

    for db_name in db_names:
        print "######################################"
        print db_name, full_name[db_name]
        switch_to_db (cursor, db_name)

        ############################
        qry = "select  hugo_symbol, chromosome, start_position, end_position, strand "
        qry += " from somatic_mutations"
        qry += " where  hugo_symbol like '%5%RNA%' "
        rows = search_db (cursor, qry)
        if not rows: continue
        for row in rows:
            print row
        print
        ############################

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

