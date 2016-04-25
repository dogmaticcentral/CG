#!/usr/bin/python

import MySQLdb
from   tcga_utils.mysql   import  *


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL",  "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",  "MESO", "OV",   "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
    
    # I don't want to start this by mistake - remove comment and put in db names as needed
    print "please comment out if you are sure you want to delete database tables"
    exit(1)
    
    for db_name in db_names:
        # check db exists
        qry = "show databases like '%s'" % db_name
        rows = search_db(cursor, qry)
        if not rows:
            print db_name, "not found"
            continue
        
        print " ** ", db_name
        switch_to_db (cursor, db_name)
        qry = "show tables"
        rows = search_db(cursor, qry)
        if not rows:
            print "no table found in", db_name
            continue
        print qry
        for row in rows:
            print "\t ", row[0]

        #for table in ( 'metastatic_mutations', 'somatic_mutations', 'mutations_meta'):
        for table in ( 'metastatic_mutations', 'somatic_mutations'):
            if ( check_table_exists (cursor, db_name, table)):
                print table, " found in ", db_name
                # if you really want to start from scratch, uncomment
                qry = "drop table %s "  % table
                rows = search_db(cursor, qry)
            else:
                print table, " not found in ", db_name
            
    cursor.close()
    db.close()
                
            

#########################################
if __name__ == '__main__':
    main()
