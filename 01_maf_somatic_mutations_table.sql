#!/usr/bin/python

# to speed up things make index:  (int 01_maf_somatic_mutations_table.sql)
# create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position)
# note: I am assuming then that no tumot will have the exact same mutation in both alleles, for the simple reason that I do nto see
# how would the thwo entried in the database then be distinguished from one another
# (rather if tumor_seq_allele1 == tumor_seq_allele2 != match_norm_seq_allele1  and tumor_seq_allele2 != match_norm_seq_allele2
# then I have somehting like that)


import MySQLdb
from   tcga_utils.mysql   import  *
import commands

 

#########################################
def make_somatic_mutations_table(cursor, db_name):

    switch_to_db (cursor, db_name)
    qry = "select database()"
    rows = search_db(cursor, qry)
    print qry
    print rows

    qry = "";
    qry += "  CREATE TABLE somatic_mutations ("
    qry += "  	 hugo_symbol VARCHAR (50) NOT NULL, "
    qry += "     entrez_gene_id INT, "
    qry += "	 aa_change BLOB, "
    qry += "	 chromosome VARCHAR (20) NOT NULL, "
    qry += "	 start_position INT  NOT NULL, "
    qry += "	 end_position INT NOT NULL, "
    qry += "	 strand VARCHAR (5) NOT NULL, "
    qry += "	 variant_classification VARCHAR (50) NOT NULL, "
    qry += "	 variant_type VARCHAR (20) NOT NULL, "
    qry += "	 reference_allele BLOB NOT NULL, "
    qry += "	 tumor_seq_allele1 BLOB NOT NULL, "
    qry += "	 tumor_seq_allele2 BLOB NOT NULL, "
    qry += "	 tumor_sample_barcode VARCHAR (50) NOT NULL, "
    qry += "	 matched_norm_sample_barcode VARCHAR (50) NOT NULL, "
    qry += "	 match_norm_seq_allele1 BLOB, "
    qry += "	 match_norm_seq_allele2 BLOB, "
    qry += "	 tumor_validation_allele1 BLOB, "
    qry += "	 tumor_validation_allele2 BLOB, "
    qry += "	 match_norm_validation_allele1 BLOB, "
    qry += "	 match_norm_validation_allele2 BLOB, "
    qry += "	 verification_status VARCHAR (20), "
    qry += "	 validation_status VARCHAR (20) NOT NULL, "
    qry += "	 mutation_status VARCHAR (50) NOT NULL"
    qry += ") ENGINE=MyISAM"
    rows = search_db(cursor, qry)
    print qry
    print rows

    qry = "";
    qry += "create index hugo_idx on somatic_mutations (hugo_symbol)"
    rows = search_db(cursor, qry)
    print qry
    print rows

    qry = "";
    qry += "create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position)" 
    rows = search_db(cursor, qry)
    print qry
    print rows



#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["COAD",  # after this we go alphabetically
                 "ACC", "BLCA", "BRCA", "CESC",  "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]


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
            make_somatic_mutations_table(cursor, db_name)

        print " ** ", db_name
        switch_to_db (cursor, db_name)
        qry = "show tables"
        rows = search_db(cursor, qry)
        print qry
        print rows
      
        table = 'somatic_mutations'

        if ( check_table_exists (cursor, db_name, table)):
            #print table, " found in ", db_name
            qry = "drop table somatic_mutations"
            rows = search_db(cursor, qry)
            make_somatic_mutations_table(cursor, db_name)
        else:
            print table, " not found in ", db_name
            exit(1)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
