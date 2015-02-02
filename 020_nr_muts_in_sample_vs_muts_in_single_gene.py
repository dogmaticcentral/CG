#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb, math
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

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS", "UVM"];

    #db_names  = ["COAD"];
   
    table = 'somatic_mutations'

    avg_size = {}
    avg_size_sq = {}
    count   = {}
    for db_name in db_names:
        print "######################################"
        print db_name, full_name[db_name]
        switch_to_db (cursor, db_name)

 
        ############################
        qry  = "select  distinct tumor_sample_barcode from somatic_mutations "
  
        rows = search_db (cursor, qry)
        if not rows: continue

        for row in rows:
            [tumor_sample_barcode] = row
            # how many mutations in this particular sample?
            qry = "select count(1) from somatic_mutations where tumor_sample_barcode = '%s'" % tumor_sample_barcode
            rows = search_db(cursor, qry)
            if not rows: # no somatic mutations (?)
                tot_number_of_mutations_in_sample = 0
                tot_number_of_mutations_in_qry   = 0
            else:
                tot_number_of_mutations_in_sample = rows[0][0]
                qry  = "select count(1) from somatic_mutations where tumor_sample_barcode = '%s'" % tumor_sample_barcode
                qry += " and hugo_symbol = '%s' " % gene_symbol
                #qry += " and not  variant_classification='silent' "
                qry += " and  variant_classification='Missense_Mutation' "
                rows2 = search_db(cursor, qry)
                if not rows2:
                    tot_number_of_mutations_in_qry   = 0
                else:
                    tot_number_of_mutations_in_qry   = rows2[0][0]

            #print tumor_sample_barcode, tot_number_of_mutations_in_sample, tot_number_of_mutations_in_qry
            if not tot_number_of_mutations_in_qry in avg_size.keys(): 
                avg_size[tot_number_of_mutations_in_qry] = 0
                avg_size_sq[tot_number_of_mutations_in_qry] = 0
                count[tot_number_of_mutations_in_qry] = 0
            
            count[tot_number_of_mutations_in_qry]    += 1
            avg_size[tot_number_of_mutations_in_qry] += tot_number_of_mutations_in_sample
            avg_size_sq[tot_number_of_mutations_in_qry] += tot_number_of_mutations_in_sample*tot_number_of_mutations_in_sample

        ############################

    for tot_number_of_mutations_in_qry, size in avg_size.iteritems():
        avg = float(size)/count[tot_number_of_mutations_in_qry]
        avg_sq = float( avg_size_sq[tot_number_of_mutations_in_qry] )/count[tot_number_of_mutations_in_qry]
        stdev = math.sqrt (avg_sq - avg*avg)
        print " %2d  %5d  %8.2f %8.2f  " % (tot_number_of_mutations_in_qry, count[tot_number_of_mutations_in_qry], avg, stdev)
    

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

