#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);


import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *

      
  
#########################################
def check_name_resolution (cursor, db_name):

    switch_to_db (cursor, db_name)
    qry = "select count(1) from name_resolution"
    rows = search_db (cursor, qry)
    if not rows or  type(rows[0][0])==str: return
    print 'total names:', rows[0][0]
    tot = float(rows[0][0])

    qry = "select count(1) from name_resolution where comment='failure'"
    rows = search_db (cursor, qry)
    if not rows or  type(rows[0][0])==str: return
        
    print '   failures:', rows[0][0]
    fail = float(rows[0][0])
    print '   fraction: %d%%' % int(fail/tot*100)
   
 
#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    db_names  = [ "ACC", "BLCA", "BRCA", "CESC", "COAD",  "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]
    for db_name  in db_names:
        print '==================================='
        print db_name
        check_name_resolution (cursor, db_name)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

