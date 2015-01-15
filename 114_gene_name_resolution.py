#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);

import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *

#########################################
def find_ensembl_id (cursor, gene):
    
    ensembl_id = ""
    comment    = ""

    if gene=='Uknown': return [ensembl_id, 'unk']

    if gene[:4] =='ENSG': 
        ensembl_id = gene
        return [ensembl_id, 'ensembl']

    # look up ens id directly
    qry  = "select distinct(ensembl_gene_id) from  baseline.hgnc_id_translation  where approved_symbol='%s'" % gene
    rows = search_db (cursor, qry)
    if rows:
        ensembl_ids = [row[0] for row in rows if row[0] > 0]
        no_hits  = len(ensembl_ids)
        if no_hits == 1:  
            ensembl_id = ensembl_ids[0].replace(' ','')
            if ensembl_id:  return [ensembl_id, 'ok']

    

    # if this didn't work, go to entrez
    qry  = "select distinct(entrez_gene_id) from  somatic_mutations  where hugo_symbol='%s'" % gene
    rows = search_db (cursor, qry)
    if rows:
        entrez_ids = [row[0] for row in rows if row[0] > 0]
        no_hits  = len(entrez_ids)
        if no_hits == 1: 
            #print gene, "single entrez id:", entrez_ids[0]
            qry  = "select distinct(ensembl_gene_id) from baseline.hgnc_id_translation  "
            qry += "where entrez_gene_id='%s'" % entrez_ids[0]
            rows = search_db (cursor, qry)
            if rows:
                ensembl_ids = [row[0] for row in rows if row[0] > 0]
                no_hits  = len(ensembl_ids)
                if no_hits == 1:   
                    ensembl_id = ensembl_ids[0].replace(' ','')
                    if ensembl_id:  return [ensembl_id, 'from hgnc_id_translation via entrez_gene_id']
                # this apparently does not happen
                #print gene, entrez_id, 'multiple ensembl', ensembl_ids

    

    # if entrez failed, try older hugo names
    for column in ['previous_symbols', 'synonyms']:
        qry = "select approved_symbol, %s, ensembl_gene_id  from baseline.hgnc_id_translation " % column
        qry += "where %s like '%%%s%%'" % (column,gene)
        rows = search_db (cursor, qry)
        if  rows:
            for row in rows:
                # check that the match is exact
                previous_symbols = row[1].replace(' ','').split(',')
                exact = [ x  for x in previous_symbols if x==gene]
                if exact:
                    hugo_id    = row[0]
                    ensembl_id = row[2].replace(' ','')
                    if ensembl_id:
                        comment =  'in hgnc_id_translation ' + column + ' for ' + hugo_id
                        return [ensembl_id, comment]

    # in the case of failure, check out ncbi translation table
    qry  = "select  symbol, synonyms, db_xrefs from  baseline.ncbi_id_translation  where  symbol='%s'" % gene
    rows = search_db (cursor, qry)
    if rows:
        db_xrefs = [row[2] for row in rows]
        ensembl_ids = []
        for db_xref in db_xrefs:
            ensembl =[ x for x in  db_xref.split('|') if 'Ensembl' in x]
            for ens in ensembl:
                ensembl_ids.append (ens.split(':')[1])
        if len(ensembl_ids)==1: 
            ensembl_id = ensembl_ids[0].replace(' ','')                    
            return [ensembl_id, "symbol in ncbi_id_translation"]

    # locus tag?
    qry  = "select  symbol, synonyms, db_xrefs from  baseline.ncbi_id_translation  where  locus_tag='%s'" % gene
    rows = search_db (cursor, qry)   
    if rows:
        db_xrefs = [row[2] for row in rows]
        ensembl_ids = []
        for db_xref in db_xrefs:
            ensembl =[ x for x in  db_xref.split('|') if 'Ensembl' in x]
            for ens in ensembl:
                ensembl_ids.append (ens.split(':')[1])
        if len(ensembl_ids)==1: 
            ensembl_id = ensembl_ids[0].replace(' ','')                    
            return [ensembl_id, "locus_tag in ncbi_id_translation"]
    

    # synonyms?
    qry  = "select  symbol, synonyms, db_xrefs from  baseline.ncbi_id_translation  where synonyms like '%%%s%%'" % gene
    rows = search_db (cursor, qry)
    if rows:
        ensembl_ids = []
        for row in rows:
            synonyms = row[1].replace(' ','').split('|')
            # check that the match is exact
            exact = [ x  for x in synonyms if x==gene]
            if exact:
                db_xrefs = [row[2] for row in rows]
                for db_xref in db_xrefs:
                    ensembl =[ x for x in  db_xref.split('|') if 'Ensembl' in x]
                    for ens in ensembl:
                        ensembl_ids.append (ens.split(':')[1])

        if len(ensembl_ids)==1: 
            ensembl_id = ensembl_ids[0]
            return [ensembl_id, "in synonyms  for " + row[0] + " in  ncbi_id_translation"]

    # try uniprot resolution table 
    fields = gene.split('.')
    fields.pop()
    gene = '.'.join(fields)
    qry  = "select * from  baseline.uniprot_id_translation  where other_db_id='%s'" % gene
    rows = search_db (cursor, qry)
    if rows:
        [uniprot_id, other_db, other_db_id ] = rows[0]
        qry  = "select * from  baseline.uniprot_id_translation  where uniprot_id='%s' " % uniprot_id
        qry += "and other_db='Enembl'" # ditto: its Enembl not Ensembl
        rows2 = search_db (cursor, qry)
        if rows2:
            [uniprot_id, other_db, other_db_id ] = rows2[0]
            ensembl_id = other_db_id
            return [ensembl_id, "uniprot_id_translation"]
      
   
    return [ensembl_id, 'failure']

      
#########################################
def make_name_resolution_table(cursor, db_name):
    qry = "";
    qry += " CREATE TABLE name_resolution ("
    qry += "  	 tcga_hugo_symbol VARCHAR (50) NOT NULL, "
    qry += "  	 ensembl_gene_id  VARCHAR (20), "
    qry += "  	 comment  BLOB "
    qry += ") ENGINE=MyISAM"
    rows = search_db(cursor, qry)
    qry = "create index hugo_idx on name_resolution (tcga_hugo_symbol)"
    rows = search_db(cursor, qry)
    print qry
    print rows
  
#########################################
def resolve_gene_name (cursor, db_name):

    switch_to_db (cursor, db_name)
   
    table = 'name_resolution'
    if ( check_table_exists (cursor, db_name, table)):
        print table, " found in ", db_name
    else:
        print 'making name_resolution table'
        make_name_resolution_table(cursor, db_name)

    table = 'somatic_mutations'
    ############################
    print "number of different genes:"
    qry = "select distinct(hugo_symbol) from somatic_mutations"
    rows = search_db(cursor, qry)
    genes = [row[0] for row in  rows]
    print "\t", len(genes)

    ############################
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
       
        [ensembl_id, comment] = find_ensembl_id (cursor, gene)
  
        fixed_fields = {'tcga_hugo_symbol':gene}
        update_fields = {'ensembl_gene_id': ensembl_id, 'comment':comment}
        if 'locus' in comment:
            store_or_update (cursor, 'name_resolution', fixed_fields, update_fields, verbose=True)
        else:
            store_or_update (cursor, 'name_resolution', fixed_fields, update_fields)
 
#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    db_names  = [  "BRCA", "ACC", "BLCA", "CESC", "COAD",  "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]
    #db_names = ['CESC']
    for db_name  in db_names:
        print '==================================='
        print db_name
        resolve_gene_name (cursor, db_name)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

