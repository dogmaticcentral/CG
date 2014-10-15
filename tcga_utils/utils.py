import MySQLdb, commands, re, os, time
from mysql import *

def hugo2ensembl (cursor, hugo_id):
    
    ensembl_id = ""
    switch_to_db(cursor, 'baseline')
    qry = "select ensembl_gene_id from hgnc_id_translation where approved_symbol='%s'" %  hugo_id
    rows  = search_db (cursor, qry)
    if rows and rows[0]:
        return rows[0][0]

    for alternative_column in ['previous_symbols','synonyms']:
        qry  = "select ensembl_gene_id, %s " %  alternative_column
        qry += "from hgnc_id_translation where ensembl_gene_id is not null and  ensembl_gene_id != ''"
        qry += "and %s like '%%%s%%'  " %  (alternative_column, hugo_id)
        rows  = search_db (cursor, qry)
    
        if rows:
            for row in rows:
                possible_ensembl  = row[0]
                alter_symbols  = row[1].replace (' ', '')
                fields = alter_symbols.split( ',')
                for field in fields:
                    if field == hugo_id:
                        return possible_ensembl
 
    return ""
