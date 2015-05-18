#!/usr/bin/python

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


    
