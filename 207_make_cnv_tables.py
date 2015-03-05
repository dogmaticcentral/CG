#!/usr/bin/python
import os
import MySQLdb
from   tcga_utils.mysql   import  *

#########################################
# watch out for  MyISAM and InnoDB tables difference!
# apparently InnoDB  is the default, which later requires db.commit() before closing the connection!

#########################################
def make_cnv_table(cursor):

    table = 'cnv_snp'
    
    qry  = "CREATE TABLE " + table + "  (cnv_id INT(10) PRIMARY KEY AUTO_INCREMENT)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

 
    for column_name in ['gene_stable_id']:
        qry = "ALTER TABLE %s  ADD %s VARCHAR(120)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False


    for column_name in ['source_ids', 'log_fold_changes', 'num_probes', 'ranges']:
        qry = "ALTER TABLE %s  ADD  %s blob " %  (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    return True


#########################################
def make_cnv_source_table(cursor):    
    table = 'cnv_source_ids'
    
    qry  = "CREATE TABLE  %s " % table
    qry += "(source_id INT(10), source_name VARCHAR(200), PRIMARY KEY (source_id) )"
    rows = search_db (cursor, qry)
    if (rows):
        print qry
        return False

    
    return True
    
    

#########################################
def main():
    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS", "UVM"]
    db     = connect_to_mysql()
    cursor = db.cursor()
    for db_name in db_names:
        print "######################################"
        print db_name
        switch_to_db (cursor,db_name)
        if not check_table_exists (cursor, db_name, 'cnv_snp'):
            print 'cnv_snp not found; making one'
            if not make_cnv_table(cursor):
                print 'error making cnv_snp'
                exit(1)
            pass
        else:
            #qry = "drop table cnv_snp"
            #search_db (cursor, qry)
            print 'cnv_snp found; creating index'
            create_index (cursor, db_name, 'stable_id_index', 'cnv_snp', ['gene_stable_id'])

        if not check_table_exists (cursor, db_name, 'cnv_source_ids'):
            print 'cnv_source_ids not found; making one'
            if not make_cnv_source_table(cursor):
                print 'error making cnv_source_ids'
                exit(1)
            pass
        else:
            #qry = "drop table cnv_source_ids"
            #search_db (cursor, qry)
            print "cnv_source_ids found"
            pass

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
