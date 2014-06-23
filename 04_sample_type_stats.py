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
    print "sample type breakdown:"
    qry  = "select TumorSampleBarcode from somatic_mutations "
    rows = search_db(cursor, qry)
    count = {}
    for  row in rows:
        tbarcode = row[0]
        # the fields are 
        # project - tissue source site (TSS)  - participant -
        # source.vial - portion.analyte  - plate - (sequenncing or charcterization center)
        fields       = tbarcode.split('-')
        sample_type  = fields[3][:2] # the first two digits fo the third field
        if not count.has_key(sample_type):  count[sample_type] = 0
        count[sample_type] += 1
        
    for sample_type, ct in count.iteritems():
        print "\t %2s   %5d " % (sample_type, ct)
    
    
            
#########################################
if __name__ == '__main__':
    main()


    
