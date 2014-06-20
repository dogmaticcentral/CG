#!/usr/bin/python

import MySQLdb
from   tcga_utils.mysql   import  *


#########################################
def check_fields(cursor, header_fields):

    for field in header_fields:
        qry  = "SELECT *  FROM information_schema.COLUMNS "
        qry += " WHERE TABLE_SCHEMA = 'COAD' AND TABLE_NAME = 'somatic_mutations'"
        qry += " and  COLUMN_NAME = '%s'" % field;
        rows = search_db (cursor, qry)
        if not rows:
            print field, " not found in somatic_mutations table", 
            exit(1)

#########################################
def store (cursor, header_fields, fields):
    
    fixed_fields  = {}
    update_fields = {}
    
    for i in range( len(header_fields) ):
        field = fields[i]
        header = header_fields[i]
        if (header in ['TumorSampleBarcode', 'Chromosome', 'Strand', 'Startposition'] ):
            fixed_fields[header] = field
        else:
            update_fields[header] = field
        
    ok = store_or_update (cursor, 'somatic_mutations', fixed_fields, update_fields)
    if not ok:
        for  i in range( len(header_fields) ):
            print  header_fields[i],  fields[i]
        exit(1)

 
#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_name  = 'COAD'
    maffile  = '/Users/ivana/databases/TCGA/COAD/Somatic_Mutations/'
    maffile += 'BCM__Mutation_Calling/Level_2/'
    maffile += 'hgsc.bcm.edu_COAD.SOLiD_DNASeq.1.somatic.maf'
    #maffile += 'hgsc.bcm.edu_COAD.IlluminaGA_DNASeq.1.somatic.maf'

    switch_to_db (cursor, db_name)

    table = 'somatic_mutations'

    if ( check_table_exists (cursor, db_name, table)):
        print table, " found in ", db_name
    else:
        print table, " not found in ", db_name
        exit(1)

    if 0: # useless bc the number of columns is not standardized
        qry  = "select count(*) from information_schema.columns "
        qry += " where table_name = 'somatic_mutations' "
        qry += " and table_schema = 'COAD' "
        rows = search_db (cursor, qry)
        number_of_columns = rows[0][0]

    inff = open(maffile, "r")
    line_ct = 0
    for line in inff:
        line = line.rstrip()
        fields = line.split ('\t')
        line_ct += 1
        if line_ct < 2: continue
        if line_ct == 2:
            header_fields = [x.replace('_', '') for x in fields]
            check_fields(cursor, header_fields)
        else:
            fields_clean = [x.replace("'", '') for x in fields]
            store (cursor, header_fields, fields_clean)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
