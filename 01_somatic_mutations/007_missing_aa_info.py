#!/usr/bin/python

import MySQLdb
from sets import Set
from   tcga_utils.mysql   import  *
import commands


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    table  = 'somatic_mutations'
    gene_names  = ['RPL5', 'RPL11','TP53']

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL",  "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",  "MESO", "OV",   "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    if False:
        for db_name in db_names:
            print " ** ", db_name
            switch_to_db (cursor, db_name)
            for gene_name in gene_names:
                qry  = "select count(1) from somatic_mutations where hugo_symbol='%s' and aa_change='missing'" %gene_name
                rows = search_db(cursor, qry)
                if not rows: continue
                for  row in rows:
                    print "\t %s"%gene_name, row[0]
        print
        print
        
    for db_name in db_names:

        switch_to_db (cursor, db_name)

        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        print db_name, table
        print "\t number of entries:", rows[0][0]
        qry  = "select count(1) from %s " % (table)
        qry += "where variant_classification='Missense_Mutation'"
        rows = search_db(cursor, qry)
        print "\t number of missense mutations:", rows[0][0]
      
        qry  = "select hugo_symbol, chromosome, start_position, end_position,  aa_change, "
        qry += " cdna_change from %s  where aa_change='missing'" % (table)
        qry += " and variant_classification='Missense_Mutation'"
        rows = search_db(cursor, qry)
        if not rows:
            number_missing = 0
        else:
            number_missing = len(rows)
        print "\t number missing the aa change info: ", number_missing
        #for  row in rows:
        #    print "\t", row
        print


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
