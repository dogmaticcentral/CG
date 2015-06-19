#!/usr/bin/python

# to speed up things make index:  (int 01_maf_somatic_mutations_table.sql)
# create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position)
# note: I am assuming then that no tumot will have the exact same mutation in both alleles, for the simple reason that I do nto see
# how would the thwo entried in the database then be distinguished from one another
# (rather if tumor_seq_allele1 == tumor_seq_allele2 != match_norm_seq_allele1  and tumor_seq_allele2 != match_norm_seq_allele2
# then I have somehting like that)


import MySQLdb
from   tcga_utils.mysql   import  *



#########################################
def make_gene_distro_table(cursor, db_name):

    table = 'rnaseq_distro_description'
    switch_to_db (cursor, db_name)
    qry = "select database()"
    rows = search_db(cursor, qry)
    print qry
    print rows

    qry = "";
    qry += "  CREATE TABLE %s ("  % table
    qry += "  	 symbol CHAR (50) NOT NULL, "
    qry += "	 number_of_points INT  DEFAULT NULL, "
    qry += "	 min FLOAT (10,2) DEFAULT NULL, "
    qry += "	 max FLOAT (10,2) DEFAULT NULL, "
    qry += "	 mean FLOAT (10,2) DEFAULT NULL, "
    qry += "	 stdev FLOAT (10,2) DEFAULT NULL, "
    qry += "	 skewness FLOAT (10,2) DEFAULT NULL, "
    qry += "	 kurtosis FLOAT (10,2) DEFAULT NULL,"
    qry += "	 distro CHAR (50) DEFAULT NULL,"
    qry += "	 KL_pval FLOAT (10,2) DEFAULT NULL, "
    qry += "	 left_cut  INT DEFAULT NULL, "
    qry += "	 right_cut INT DEFAULT NULL, "
    qry += "     shape FLOAT (10,2) DEFAULT NULL,"
    qry += "     location FLOAT (10,2) DEFAULT NULL,"
    qry += "     scale FLOAT (10,2) DEFAULT NULL,"
    qry += "     interval_endpoints BLOB"
    qry += ") ENGINE=MyISAM"
    rows = search_db(cursor, qry)
    print qry
    print rows

    if True:
        qry = "";
        qry += "create index hugo_idx on %s (symbol)" % table
        rows = search_db(cursor, qry)
        print qry
        print rows



#########################################
def modify_gene_distro_table(cursor, db_name):

    switch_to_db (cursor, db_name)
    qry = "select database()"
    rows = search_db(cursor, qry)
    print qry
    print rows

    #qry = "alter table gene_expression add  source  VARCHAR (50) DEFAULT NULL"
    #rows = search_db(cursor, qry)
    #print qry
    #print rows
    
    qry = "update gene_expression set source = 'unc' "
    rows = search_db(cursor, qry)
    print qry
    print rows
    


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["BLCA","BRCA","COAD", "HNSC", "KIRC","KIRP","LIHC","LUAD","LUSC","REA","UCEC"]

    for db_name in db_names:
        # check db exists
        qry = "show databases like '%s'" % db_name
        rows = search_db(cursor, qry)
        if not rows:
            print db_name, "not found"
            qry = "create database " +  db_name
            rows = search_db(cursor, qry)
            print qry
            print rows

        print " ** ", db_name
        switch_to_db (cursor, db_name)
        qry = "show tables"
        rows = search_db(cursor, qry)
        print qry
        print rows

        table = 'rnaseq_distro_description'

        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
            # if you really want to start from scratch, uncomment
            #qry = "drop table %s" % table
            #rows = search_db(cursor, qry)
            #make_gene_distro_table(cursor, db_name)
            #modify_gene_distro_table(cursor, db_name)
        else:
            print table, " not found in ", db_name
            make_gene_distro_table(cursor, db_name)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
