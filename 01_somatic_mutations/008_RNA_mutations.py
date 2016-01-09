#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *



#########################################
def main():


    full_name = read_cancer_names ()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG",  "LIHC", "LUAD", "LUSC", "OV",  "PAAD", "PCPG", "PRAD", "REA",
                 "SARC","SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    db_names = ["LUAD"]

   
    table = 'somatic_mutations'

    for db_name in db_names:
        print "######################################"
        print db_name, full_name[db_name]
        switch_to_db (cursor, db_name)

        ############################
        qry = "select  hugo_symbol, chromosome, start_position, end_position, strand "
        qry += " from somatic_mutations"
        qry += " where  hugo_symbol like '%5%RNA%' "
        rows = search_db (cursor, qry)
        if not rows: continue
        for row in rows:
            print row
        print
        ############################

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

