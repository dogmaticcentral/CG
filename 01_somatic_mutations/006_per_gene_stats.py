#!/usr/bin/python -u
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *
from   tcga_utils.ensembl   import  *

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    genes = ["AKAP13", "ESR1", "HAND2", "PRKACA", "PRKAR2A", "PRKAR2B", "PRKCA"]
    full_name = read_cancer_names ()
    table = 'somatic_mutations'
    for db_name in db_names:
        print "######################################"
        print db_name, full_name[db_name]
        switch_to_db (cursor, db_name)


        ############################
        print "number of entries:",
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        print rows[0][0]
        ############################

        ############################
        print "number of patients:",
        qry = "select distinct(sample_barcode_short) from somatic_mutations"
        rows = search_db(cursor, qry)
        patients = [row[0] for row in  rows]
        total_patients = len(patients)
        print  total_patients


        ############################
        if not genes: # go for all of them
            print "number of different genes:"
            qry = "select distinct(hugo_symbol) from somatic_mutations"
            rows = search_db(cursor, qry)
            genes = [row[0] for row in  rows]
            print "\t", len(genes)

        ############################
        print "mutations reported per gene"
        print " %10s   %5s  %5s    %s " % ("gene_name", "silent", "non_silent", "silent/non")
        for gene in genes:

            [silent_ct, non_silent_ct] = silent_proportion(cursor, gene)
            if non_silent_ct:
                 print " %10s   %5d  %5d    %4.2f  " % (gene, silent_ct,
                                                                 non_silent_ct, float(silent_ct)/non_silent_ct)
            else:
                print " %10s   %5d  %5d  all_silent " % (gene, silent_ct, non_silent_ct)

            #print "  %4d  %10s   %5d  %5d  " % ( ct, gene, entries_per_gene[gene],   silent_per_gene[gene])
        print


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

