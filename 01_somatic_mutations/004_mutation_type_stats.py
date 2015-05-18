#!/usr/bin/python

import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *

#########################################
def main():

    extra_genes = ["RPL5", "RPL11"]
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    table = 'somatic_mutations'

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


    full_name = read_cancer_names ()


    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    # the order in which we want the variants output:
    variant_order = ["Missense_Mutation", "Silent", "Nonsense_Mutation", "RNA", "Splice_Site", "Frame_Shift_Ins", 
                     "Frame_Shift_Del", "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation", 
                     "5UTR", "3UTR", "IGR", "Intron", "5Flank"]
    grand_total = {}
    for variant in variant_order:
        grand_total[variant] = 0

    grand_total_extra = {}
    for gene in extra_genes: 
        grand_total_extra[gene] = {}
        for variant in variant_order:
            grand_total_extra[gene][variant] = 0
        

    total_extra    = {}
    total_patients = 0
    for db_name in db_names:

        total = 0
        for gene in extra_genes: 
            total_extra[gene] = 0
        ############################
        switch_to_db (cursor, db_name)
        print 
        print "################################"
        print db_name, full_name[db_name]
  
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        total = int (rows[0][0])
        print "number of entries:", total

        qry = "select distinct sample_barcode_short from somatic_mutations"
        rows = search_db(cursor, qry)
        number_of_patients = len(rows)
        total_patients    += number_of_patients
        print "number of patients:", number_of_patients
        ############################

        for gene in extra_genes: 
            qry = "select count(1) from " + table
            qry += " where hugo_symbol = '%s' " % gene
            rows = search_db(cursor, qry)
            total_extra[gene] = int (rows[0][0])
            print "\tnumber of entries for", gene, ":",  total_extra[gene]
            


        ############################
        print "variant classification (# cases):  %20s" %  "overall",
        for gene in extra_genes:
            print " %20s " % gene,
        print
        qry = "select distinct(variant_classification) from somatic_mutations"
        rows = search_db(cursor, qry)
        variants = [row[0] for row in  rows]
        for variant in variant_order:
            if not variant in variants: continue
            qry = "select count(1) from somatic_mutations where variant_classification='%s'" % variant
            rows = search_db(cursor, qry)
            print "\t %30s    %6d   (%4.1f%%)" %  (variant, rows[0][0], float(rows[0][0])/total*100 ),
            grand_total[variant] += rows[0][0]
            for gene in extra_genes:
                if not total_extra[gene]:
                    print "\t %6d   (%4.1f%%)" %  (0, 0.0),
                else:
                    qry = "select count(1) from somatic_mutations "
                    qry += " where hugo_symbol = '%s' " % gene
                    qry += " and variant_classification='%s'" % variant
                    rows = search_db(cursor, qry)
                    print "\t %6d   (%4.1f%%)" % (rows[0][0], float(rows[0][0])/total_extra[gene]*100 ),
                    grand_total_extra[gene][variant] += rows[0][0]
            print

    grand_grand_total = 0
    for variant in variant_order:
        grand_grand_total += grand_total[variant] 

    grand_grand_total_extra = {}
    for gene in extra_genes:
        grand_grand_total_extra[gene] = 0
        for variant in variant_order:
            grand_grand_total_extra[gene] += grand_total_extra[gene][variant] 
       

    print 
    print "################################"
    print 'pan-cancer'
    print "number of entries:", grand_grand_total
    print "number of patients:", total_patients
    print "variant classification (# cases):  %20s" %  "overall",
    for gene in extra_genes:
        print " %20s " %  gene,
    print
    for variant in variant_order:
        print "\t %30s    %8d   (%4.1f%%)" %  (variant,  grand_total[variant], float(grand_total[variant] )/grand_grand_total*100 ),
        for gene in extra_genes:
            if not grand_grand_total_extra[gene]:
                print "\t %6d   (%4.1f%%)" %  (0, 0.0),
            else:
                qry = "select count(1) from somatic_mutations "
                print "\t %6d   (%4.1f%%)" %  (grand_total_extra[gene][variant], 
                                               float(  grand_total_extra[gene][variant] )/grand_grand_total_extra[gene]*100 ),
        print
        
  

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

