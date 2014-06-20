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

    switch_to_db (cursor, db_name)

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    ############################
    print db_name, table, "number of entries:"
    qry = "select count(1) from " + table
    rows = search_db(cursor, qry)
    print "\t", rows[0][0]
    print 
    ############################

    ############################
    print "variant classification:"
    qry = "select distinct(variantclassification) from somatic_mutations"
    rows = search_db(cursor, qry)
    variants = [row[0] for row in  rows]
    for variant in variants:
        qry = "select count(1) from somatic_mutations where variantclassification='%s'" % variant
        rows = search_db(cursor, qry)
        print "\t", variant, rows[0][0], "cases"

 

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

