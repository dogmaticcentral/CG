#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *

#########################################
def find_ensembl_id (cursor, gene):
    
    ensembl_id = ""

    qry  = "select distinct(ensembl_gene_id) from  baseline.id_translation  where approved_symbol='%s'" % gene
    rows = search_db (cursor, qry)
    if rows:
        ensembl_ids = [row[0] for row in rows if row[0] > 0]
        no_hits  = len(ensembl_ids)
        if no_hits == 1:   return ensembl_ids[0]
  
 
    # if this didn't work, go to entrez
    qry  = "select distinct(entrez_gene_id) from  baseline.id_translation  where approved_symbol='%s'" % gene
    if rows:
        rows = search_db (cursor, qry)
        entrez_ids = [row[0] for row in rows if row[0] > 0]
        no_hits  = len(entrez_ids)
        if no_hits == 1: 
            #print gene, "single entrez id:", entrez_ids[0]
            qry  = "select distinct(ensembl_id) from baseline.id_translation  "
            qry += "where entrez_gene_id='%s'" % entrez_ids[0]
            rows = search_db (cursor, qry)
            ensembl_ids = [row[0] for row in rows if row[0] > 0]
            no_hits  = len(ensembl_ids)
            if no_hits == 1:   return ensembl_ids[0]
            # this apparently does not happen
            print gene, entrez_id, 'multiple ensembl', ensembl_ids

    # if entrez failed, try aolder hugo names
    for column in ['previous_symbols', 'synonyms']:
        qry = "select approved_symbol, %s, ensembl_gene_id  from baseline.id_translation " % column
        qry += "where %s like '%%%s%%'" % (column,gene)
        rows = search_db (cursor, qry)
        if  rows:
            for row in rows:
                # check that the match is exact
                previous_symbols = row[1].replace(' ','').split(',')
                exact = [ x  for x in previous_symbols if x==gene]
                if exact:
                    hugo_id    = row[0]
                    ensembl_id =  row[2]
                    if ensembl_id:
                        #print gene, 'found in', column, 'for', hugo_id, 'ensembl id', row[2]
                        return ensembl_id
    
    print 'not found', gene
    return ensembl_id

      
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
    qry = "select distinct(hugo_symbol) from somatic_mutations"
    rows = search_db(cursor, qry)
    genes = [row[0] for row in  rows]
    print "\t", len(genes)

    ############################
    print "mutations reported per gene"
    entries_per_gene = {}
    silent_per_gene  = {}
    for gene in genes:
        switch_to_db (cursor, db_name)
        qry = "select count(1) from somatic_mutations where hugo_symbol='%s'" % gene
        rows = search_db(cursor, qry)
        entries_per_gene[gene] = rows[0][0]
        qry  += " and variant_classification='silent' "        
        rows = search_db(cursor, qry)
        silent_per_gene[gene] = rows[0][0]
        # what is the length of the peptide this gene translates to?
       
        ensembl_id = find_ensembl_id (cursor, gene)
        if ensembl_id:
            #print "ensembl id found for ", ensembl_id, gene
            pass
        continue
        qry = "select peptide from baseline.aa_sequence where  ensembl_id='%s'" % gene
        rows = search_db(cursor, qry)
        if not rows:
            print 'no peptide found for ', gene
            exit(1)
        # total peptide length
    exit(1)

    sorted_genes =  sorted(genes, key= lambda x: -entries_per_gene[x])
    ct = 0
    for gene in sorted_genes:
        ratio = float(silent_per_gene[gene])/entries_per_gene[gene]
        #if ratio > 0.25 or entries_per_gene[gene] < 5: continue
        ct += 1
        #print " %4d   %10s   %5d  %5d    %4.2f   %5.2f%%" % (ct, gene, entries_per_gene[gene],
        #                                                     silent_per_gene[gene], ratio,  float(ct)/len(genes)*100)
        print "  %4d  %10s   %5d   %8.2f  %5d  " % ( ct, gene, entries_per_gene[gene],  expected[gene], silent_per_gene[gene])
    print


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

