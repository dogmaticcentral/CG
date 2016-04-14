import MySQLdb, commands, re, os, time
from mysql import *

################################################################################################
def process_header_line(headerline):
    
    if len(headerline)==0: return []

    header_fields = [x.lower() for x in headerline.split('\t')]
    # translate the header nomenclature to what we expect to see:
    for i in range(len(header_fields)):
        if header_fields[i] == 'chrom':
            header_fields[i] ='chromosome'
        elif header_fields[i] in ['amino_acid_change_wu',
                              'aachange', 'amino_acid_change',
                              'protein_change']:
            header_fields[i] = 'aa_change'
        elif header_fields[i] in ['cdna_change', 'chromchange', 'c_position_wu', 'c_position']:
            header_fields[i] = 'cdna_change'

    # I am adding this one so I do not have to search the database by doing substring comparison
    header_fields.append('sample_barcode_short')

    return header_fields

################################################################################################
# Required_fields are the absolute minimum we need to 
# reconstruct the mutation - that they are missing  should not happen at all
def get_required_fields ():
    return ['start_position',  'tumor_seq_allele1', 'tumor_seq_allele2',
            'match_norm_seq_allele1',  'match_norm_seq_allele2']
 
################################################################################################
# The fields we have predeclared in the database are what we expect to see in the file
# The rest of the code should handle the missing fields somehow.
def get_expected_fields(cursor, db_name, table):

    qry  = "describe %s " % table
    rows = search_db (cursor, qry)
    if not rows:
        print "%s not found in  %s" % (table, db_name)
        exit(1)
    # the first header field is 'id' - this is entry id in the table
    expected_fields = [row[0].lower() for row in rows[1:]]
    return expected_fields

################################################################################################
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


#########################################
def read_cancer_names ():
    full_name= {}
    inf = open ("/Users/ivana/pypeworks/tcga/cancer_names.txt", "r")
    for line in inf:
        line   = line.rstrip()
        field = line.split ("\t")
        if field[0] == 'READ':
            field[0] = 'REA'
        full_name[field[0]] = field[1]
    inf.close()

    return full_name

