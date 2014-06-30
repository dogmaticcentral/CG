#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *

#########################################
def read_cancer_names ():
    full_name= {}
    inf = open ("/Users/ivana/pypeworks/tcga/cancer_names.txt", "r")
    for line in inf:
        line   = line.rstrip() 
        field = line.split ("\t")
        if field[0] == 'READ':
            field[0] = 'REA'
        full_name[field[0]] = field[1]
    inf.close()

    return full_name
    

#########################################
def main():

    if len(sys.argv) < 2:
        print  "usage: %s <gene symbol> " % sys.argv[0]
        exit(1)

    full_name = read_cancer_names ()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    gene_symbol = sys.argv[1].upper()
    print gene_symbol
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "COAD", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]
   
    table = 'somatic_mutations'

    for db_name in db_names:
        print "######################################"
        print db_name, full_name[db_name]
        switch_to_db (cursor, db_name)

        ############################
        print "\t number of entries:", 
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        print  rows[0][0]

        if not rows[0][0]: continue

        ############################
        print "\t mutations found in", gene_symbol + ":"
        qry = "select  hugo_symbol, variant_classification, aa_change, tumor_sample_barcode"
        qry += " from somatic_mutations"
        qry += " where  hugo_symbol = '%s' " % gene_symbol
        #qry += " and not  variant_classification='silent' "
        rows = search_db (cursor, qry)
        if not rows: continue
        for row in rows:
            print "\t\t", row
        print 
        ############################

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

