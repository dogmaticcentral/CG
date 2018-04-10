#!/usr/bin/python

#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2016 Ivana Mihalek.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#

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
def make_gene_expression_table(cursor, db_name):

    switch_to_db (cursor, db_name)
    qry = "select database()"
    rows = search_db(cursor, qry)
    print qry
    print rows

    qry = "";
    qry += "  CREATE TABLE gene_expression ("
    qry += "  	 symbol VARCHAR (50) NOT NULL, "
    qry += "	 sample_id VARCHAR (200) DEFAULT NULL, "
    qry += "	 fold_change FLOAT (10,2) DEFAULT NULL, "
    qry += "	 source VARCHAR (50) DEFAULT NULL"
    qry += ") ENGINE=MyISAM"
    rows = search_db(cursor, qry)
    print qry
    print rows

    qry = "";
    qry += "create index hugo_idx on gene_expression (symbol)"
    rows = search_db(cursor, qry)
    print qry
    print rows
    qry = "";
    qry += "create index tumor_id_idx on gene_expression (symbol, sample_id)"
    rows = search_db(cursor, qry)
    print qry
    print rows



#########################################
def modify_gene_expression_table(cursor, db_name):

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
  
        print " ** ", db_name
        switch_to_db (cursor, db_name)
        qry = "show tables"
        rows = search_db(cursor, qry)
        print qry
        print rows
      
        table = 'gene_expression'

        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
            # if you really want to start from scratch, uncomment
            # qry = "drop table gene_expression"
            # rows = search_db(cursor, qry)
            #make_gene_expression_table(cursor, db_name)
            modify_gene_expression_table(cursor, db_name)
        else:
            print table, " not found in ", db_name
            exit(1)
            #make_gene_expression_table(cursor, db_name)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
