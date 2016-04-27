#!/usr/bin/python

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

import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_name  = 'COAD'
    table = 'somatic_mutations'

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA", "FPPP", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    for db_name in db_names:
        print " ** ", db_name
        switch_to_db (cursor, db_name)

        ############################
        print "sample type breakdown:"
        qry  = "select distinct sample_barcode_short from somatic_mutations "
        rows = search_db(cursor, qry)
        count = {}
        for  row in rows:
            short_barcode = row[0]
            source  = short_barcode[-2:] # the first two digits fo the third field
            if not count.has_key(source):  count[source] = 0
            count[source] += 1

        for source, ct in count.iteritems():
            print "\t %2s   %5d " % (source, ct)
    
    cursor.close()
    db.close()


            
#########################################
if __name__ == '__main__':
    main()


    
