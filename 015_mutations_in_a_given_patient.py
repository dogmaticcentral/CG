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

    if len(sys.argv) < 3:
        print  "usage: %s <db_name> <patient_code> " % sys.argv[0]
        exit(1)

    full_name = read_cancer_names ()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    db_name = sys.argv[1].upper()
    patient_code = sys.argv[2].upper()
    print db_name, patient_code
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    table = 'somatic_mutations'

    print "######################################"
    print db_name, full_name[db_name]
    switch_to_db (cursor, db_name)

    ############################
    qry = "select  hugo_symbol, variant_classification, aa_change, tumor_sample_barcode"
    qry += " from somatic_mutations"
    qry += " where  tumor_sample_barcode like '%%%s%%' " % patient_code
    #qry += " and not  variant_classification='silent' "
    rows = search_db (cursor, qry)
    if not rows[0][0]: 
        print qry
        print 'no return'
        exit(1)

    ############################
    print "\t number of entries for", patient_code+":", len(rows)
   

    if not rows[0][0]: 
        print '\t no entries found'
        exit(1)


    for row in rows:
        print "\t\t", row
    print 
    ############################

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

