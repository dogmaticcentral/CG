#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *



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

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]


   
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
            [hugo_symbol, variant_classification, aa_change, tumor_sample_barcode] = row
            # how many mutations in this particular sample?
            qry = "select count(1) from somatic_mutations where tumor_sample_barcode = '%s'" % tumor_sample_barcode
            rows = search_db(cursor, qry)
            if not rows: 
                tot_number_of_mutations_in_sample = 0
            else:
                tot_number_of_mutations_in_sample = rows[0][0]
            #if aa_change and 'P72' in aa_change:
            print "\t\t %10s  %20s  %10s  %6d" % ( hugo_symbol, variant_classification, aa_change, 
                                                   tot_number_of_mutations_in_sample)
        print 
        ############################

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

