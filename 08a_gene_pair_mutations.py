#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);

# extracting mutations ('catch' is the output from this script)
# grep Missense catch | grep -v found | grep RPL5 | grep -v Silent | awk '{print $5}' | sed 's/p\.//g' | sed 's/[A-Z]//g' | awk '{printf "%d+", $1}' && echo


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

    if len(sys.argv) < 3:
        print  "usage: %s <gene symbol 1>  <gene symbol 2>" % sys.argv[0]
        exit(1)

    full_name = read_cancer_names ()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    gene_symbol_1 = sys.argv[1].upper()
    gene_symbol_2 = sys.argv[2].upper()
    
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
        print "total number of entries:", 
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        print  rows[0][0]

        if not rows[0][0]: continue

        ############################
        qry = "select  variant_classification, aa_change, tumor_sample_barcode"
        qry += " from somatic_mutations"
        qry += " where  hugo_symbol = '%s' " % gene_symbol_1
        #qry += " and not  variant_classification='silent' "
        rows = search_db (cursor, qry)
        if not rows: 
            print 'no mutations found in ', gene_symbol_1
            continue
        print "%15s  %10s  %20s  %15s" % ('#muts_in_sample', 'name1', 'variant1', 'aa_change1'),
        print "%10s  %20s  %15s "  % ('name2', 'variant2', 'aa_change2')
        for row in rows:
            [variant_classification, aa_change, tumor_sample_barcode] = row
            # how many mutations in this particular sample?
            qry = "select count(1) from somatic_mutations where tumor_sample_barcode = '%s'" % tumor_sample_barcode
            rows = search_db(cursor, qry)
            if not rows: 
                tot_number_of_mutations_in_sample = 0
            else:
                tot_number_of_mutations_in_sample = rows[0][0]
            print "%s  %15d  %10s  %20s  %15s  " % (tumor_sample_barcode, tot_number_of_mutations_in_sample, gene_symbol_1,
                                               variant_classification, aa_change),
            # do we have the second gene mutated in the same sample?
            qry = "select  variant_classification, aa_change"
            qry += " from somatic_mutations"
            qry += " where  hugo_symbol = '%s' " % gene_symbol_2
            qry += " and  tumor_sample_barcode  = '%s' " % tumor_sample_barcode
            rows = search_db (cursor, qry)
            if not rows: 
                print "   %s  %20s " %  (gene_symbol_2, 'No mutations found')
            else:
                [variant_classification_2, aa_change_2] = rows[0]
                print "   %s  %20s  %15s" % ( gene_symbol_2, variant_classification_2, aa_change_2)
                for row in rows[1:]:
                    [variant_classification_2, aa_change_2] = row
                    print " "*70,
                    print " %s  %20s  %15s" % ( gene_symbol_2, variant_classification_2, aa_change_2)
        print 
        ############################

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

