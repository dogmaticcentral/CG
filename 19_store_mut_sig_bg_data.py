#!/usr/bin/python


import MySQLdb
from   tcga_utils.mysql   import  *

#########################################
def make_table(cursor, db_name, table):
    # this is three functions in one:
    if table == 'cpg_nucleotides':
        qry = "CREATE TABLE %s  (" % table
        qry += "ensembl_id  VARCHAR(20) NOT NULL KEY, "
        qry += "cpgs  LONGBLOB"

    elif table == 'possible_mutations':
        qry  = "CREATE TABLE %s  (" % table
        qry += "ensembl_id  VARCHAR(20) NOT NULL, "
        qry += "category VARCHAR(8), "
        qry += "silent    INT, "
        qry += "nonsense  INT, "
        qry += "missense  INT "
        
    elif table == 'canonical_sequence':
        qry = "CREATE TABLE %s  (" % table
        qry += "ensembl_id  VARCHAR(20) NOT NULL KEY, "
        qry += "sequence LONGBLOB"

    else:
        print 'table name', table, 'was not expected here'
        exit(1)


    qry += ") ENGINE=MyISAM"

    rows = search_db(cursor, qry)
    print qry
    print rows

    qry = "";
    qry += "create index ensembl_idx on %s  (ensembl_id)" % table
    rows = search_db(cursor, qry)
    print qry
    print rows



#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    input_file = '/Users/ivana/databases/mut_significance_bg_data.txt'

    db_name = 'baseline'
    switch_to_db (cursor, db_name)


    # check that we have all tables we need
    # otherwise make them
    for table in ('cpg_nucleotides', 'possible_mutations', 'canonical_sequence'):
        if not check_table_exists (cursor, db_name, table):
            print table, 'not found'
            make_table(cursor, db_name, table)
        else:
            print table, 'found'
            #qry = 'describe  %s' % table
            #rows = search_db(cursor, qry)
            #for row in rows:
            #    print row
           

    inf = open(input_file, "r")
    ensembl_id = ""

    [readingCpg, readingPossibleMutations, readingCanonicalSequence] = [False, False, False]
    for line in inf:
        if  line[:4] == "ENSG":
            if  'done' in line:
                ensembl_id = ""
                pass
            else:
                fields = line.split(' ')
                ensembl_id = fields[0]
                [readingCpg, readingPossibleMutations, readingCanonicalSequence] = [False, False, False]

        elif line[0] == '#':
            if line[2:5]=='CpG':
                [readingCpg, readingPossibleMutations, readingCanonicalSequence] = [True, False, False]
            elif line [2:11] == 'mutations':
                [readingCpg, readingPossibleMutations, readingCanonicalSequence] = [False, True, False]
            elif line [2:11] == 'canonical':
                [readingCpg, readingPossibleMutations, readingCanonicalSequence] = [False, False, True]
               
        else:
            fields = line.rstrip().split()
            if readingCpg:
                if  len(fields) != 1: continue
                [cpgs] = fields
                # store or update <<< HERE
                fixed_fields  = {'ensembl_id':ensembl_id}
                update_fields = {'cpgs':cpgs}
                store_or_update (cursor, 'cpg_nucleotides', fixed_fields, update_fields)

            elif readingPossibleMutations:
                if  len(fields) != 4: continue
                [category, silent, nonsense,  missense] = fields
                fixed_fields  = {'ensembl_id':ensembl_id, 'category':category}
                update_fields = {'silent':silent, 'nonsense':nonsense,  'missense':missense}
                store_or_update (cursor, 'possible_mutations', fixed_fields, update_fields)

            elif readingCanonicalSequence:
                if  len(fields) != 1: continue
                [sequence] = fields 
                fixed_fields  = {'ensembl_id':ensembl_id}
                update_fields = {'sequence':sequence}
                store_or_update (cursor, 'canonical_sequence', fixed_fields, update_fields)
            
    inf.close()
    

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
