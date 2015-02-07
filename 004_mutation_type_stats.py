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

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS", "UVM"]

    # the order in which we want the variants output:
    variant_order = ["Missense_Mutation", "Silent", "Nonsense_Mutation", "RNA", "Splice_Site", "Frame_Shift_Ins", 
                     "Frame_Shift_Del", "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation", 
                     "5UTR", "3UTR", "IGR", "Intron", "5Flank"]
    
    for db_name in db_names:
        ############################
        switch_to_db (cursor, db_name)
        print 
        print "################################"
        print db_name, table
  
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        total = int (rows[0][0])
        print "number of entries:", total
        ############################



        ############################
        print "variant classification:"
        qry = "select distinct(variant_classification) from somatic_mutations"
        rows = search_db(cursor, qry)
        variants = [row[0] for row in  rows]
        for variant in variant_order:
            if not variant in variants: continue
            qry = "select count(1) from somatic_mutations where variant_classification='%s'" % variant
            rows = search_db(cursor, qry)
            print "\t %30s    %6d cases   (%4.1f%%)" %  (variant, rows[0][0], float(rows[0][0])/total*100 )


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

