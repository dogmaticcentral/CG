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
def find_existing_fields(cursor, db_name):

    qry  = "describe somatic_mutations"
    rows = search_db (cursor, qry)
    if not rows:
        print "somatic_mutations not found in ", db_name
        exit(1)
    existing_fields = [row[0] for row in rows]
    return existing_fields 

#########################################
def store (cursor, header_fields, fields, usable_field):
    
    fixed_fields  = {}
    update_fields = {}
    
    for i in range( len(header_fields) ):
        if not usable_field[i]: continue
        field = fields[i]
        header = header_fields[i]
        if (header in ['tumor_sample_barcode', 'chromosome', 'strand', 'start_position'] ):
            fixed_fields[header] = field
        else:
            update_fields[header] = field
        
    ok = store_or_update (cursor, 'somatic_mutations', fixed_fields, update_fields)
    if not ok:
        print 'store failure:'
        print fixed_fields[header]
        print update_fields[header] 
        exit(1)

 
#########################################
def load_maf (cursor, db_name, existing_header_fields, maffile):
    inff = open(maffile, "r")
    line_ct = 0
    for line in inff:
        if line[0]=='#': continue
        line = line.rstrip()
        fields = line.split ('\t')
        line_ct += 1
        if line_ct == 1:
            
            header_fields = [x.lower() for x in fields]
            for i in range(len(header_fields)):
                if header_fields[i] == 'chrom':
                    header_fields[i] ='chromosome'
                elif header_fields[i] in ['amino_acid_change_wu',
				     'aachange', 'amino_acid_change',
				     'protein_change']:
                    header_fields[i] = 'aa_change'
            
            usable_field = map (lambda x: x in existing_header_fields, header_fields)
            
        else:
            fields_clean = [x.replace("'", '') for x in fields]
            store (cursor, header_fields, fields_clean, usable_field)
    inff.close()

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

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "COAD", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]

    # there seems to have been a problem with OV ... re-run

    db_names = ["OV"]

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
      
        db_dir  = '/Users/ivana/databases/TCGA/'+db_name
        ret      = commands.getoutput('find '+db_dir+' -name "*.maf"')
        maf_files = ret.split('\n')

        switch_to_db (cursor, db_name)

        table = 'somatic_mutations'

        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name
            exit(1)

        existing_header_fields = find_existing_fields(cursor, db_name)

        for maffile in maf_files:
            print '\t loading:', maffile
            load_maf (cursor, db_name, existing_header_fields, maffile)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
