#!/usr/bin/python

import MySQLdb
from   tcga_utils.mysql   import  *
import commands


#########################################



#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_name = "baseline"
    switch_to_db (cursor, db_name)

    # hugo2ensembl
    infile = "/Users/ivana/databases/hgnc2ensembl.txt"
    inf    = open (infile, "r")
    for line in inf:
        if 'Approved' in line: continue
        line = line.rstrip()
        fields = line.split ("\t");
        if len(fields) < 2 or not 'ENSG' in fields[1]: continue
        fixed_field  = {'hugo_symbol':fields[0].replace (' ','')}
        update_field = {'ensembl_id':fields[1].replace (' ','')}
        ok = store_or_update (cursor, 'aa_sequence', fixed_field, update_field, verbose = False)
        if not ok:
            print 'store failure:'
            print fixed_field
            print update_field 
            exit(1)
    inf.close()

    cursor.close()
    db.close()

#########################################
if __name__ == '__main__':
    main()
