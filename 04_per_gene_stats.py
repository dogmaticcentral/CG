#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_name  = 'COAD'
    table = 'somatic_mutations'

    switch_to_db (cursor, db_name)

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    ############################
    print db_name, table, "number of entries:"
    qry = "select count(1) from " + table
    rows = search_db(cursor, qry)
    print "\t", rows[0][0]
    print 
    ############################

    ############################
    print "number of different genes:"
    qry = "select distinct(hugoSymbol) from somatic_mutations"
    rows = search_db(cursor, qry)
    genes = [row[0] for row in  rows]
    print "\t", len(genes)

    ############################
    print "mutations reported per gene"
    entries_per_gene = {}
    silent_per_gene  = {}
    for gene in genes:
        qry = "select count(1) from somatic_mutations where hugoSymbol='%s'" % gene
        rows = search_db(cursor, qry)
        entries_per_gene[gene] = rows[0][0]
        qry  += " and variantclassification='silent' "        
        rows = search_db(cursor, qry)
        silent_per_gene[gene] = rows[0][0]

    sorted_genes =  sorted(genes, key= lambda x: entries_per_gene[x])
    for gene in sorted_genes:
        ratio = float(silent_per_gene[gene])/entries_per_gene[gene]
        if ratio > 0.15: continue
        print " %10s   %5d  %5d    %4.2f " % (gene, entries_per_gene[gene], silent_per_gene[gene], ratio)
    print


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

