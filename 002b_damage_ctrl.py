#!/usr/bin/python


# damage control: 
# The test for exisiting entries below
#        if (header in ['tumor_sample_barcode', 'chromosome', 'strand', 'start_position'] ):
#            fixed_fields[header] = field
# may fail, because somesubmitters may call the strand '+', and some may call it '+1'
# the sumbols used to denote strand:  ['+', '+1', '-1', '1', '-', '']
# some put in empty string, which might be right: isf a mutation exists, it will be present in both
# strands - should not have used strand as the criterion for uniqueness at all

# but then some good came out of it, bcs as it turns out, some centers deposit the info about the involved AA mutation,
# which makes my life significantly easier

# There is a serious issue that I am shoving under the carpet here:
# what if the two sets of somatic mutations are significantly different?
# It shouldn't be my call which one is different. Therefore I should probably throw away the data.
# I am solving it by pretending that the problem does not exist: the wo sets complement each other, and the difference
# is negligible.

# to speed up things make index:  (int 01_maf_somatic_mutations_table.sql)
# create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position)
# note: I am assuming then that no tumot will have the exact same mutation in both alleles, for the simple reason that I do nto see
# how would the thwo entried in the database then be distinguished from one another
# (rather if tumor_seq_allele1 == tumor_seq_allele2 != match_norm_seq_allele1  and tumor_seq_allele2 != match_norm_seq_allele2
# then I have somehting like that)




import MySQLdb
from sets import Set
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
def store (cursor, header_fields, fields):
    
    fixed_fields  = {}
    update_fields = {}

    # for some fiels, the length of the header row and the 
    # data rows is not equal
    # there isn't much I can do about it:
    shorter =  len(header_fields) 
    if  len(fields) <  shorter: 
        shorter =  len(fields) 
    
    for i in range(shorter ):
        field = fields[i]
        header = header_fields[i]
        if (header in ['tumor_sample_barcode', 'chromosome',  'start_position'] ):
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
def check_strand_annotation (cursor, db_names):
    annot_set = []
    for db_name in db_names:
        print " ** ", db_name
        switch_to_db (cursor, db_name)
        qry = "select distinct strand from somatic_mutations"
        rows = search_db(cursor, qry)
        
        for row in rows:
            for strand_annot in row:
                print strand_annot
                if strand_annot not in annot_set:
                    annot_set.append (strand_annot)
    print "annot used:", annot_set

    
#########################################
def cleanup (cursor, header_fields, tbarcode, chrom, start):
    
    
    qry = "select * from somatic_mutations where "
    qry += "tumor_sample_barcode = '%s' " % tbarcode
    qry += "and chromosome = '%s'  " % chrom
    qry += "and start_position = %s  " % start
    rows = search_db(cursor, qry)
    
    new_row   = []
    first_row = rows[0]

    no_cols = len(header_fields)
    for row in rows[1:]:
        if len (row) != no_cols:
            print "oink ?!"
            exit(1)

    for i in range(no_cols):
        field_set = False
        for row in rows:
            if str(row[i]).lower() != 'missing':
                new_row.append(row[i])
                field_set = True
                break
        if not field_set:
            new_row.append('missing')

    for i in range(no_cols):
        print i, header_fields[i], 
        for row in rows:
            print row[i],
        print new_row[i]
    # delete old rows
    qry = "delete from somatic_mutations where "
    qry += "tumor_sample_barcode = '%s' " % tbarcode
    qry += "and chromosome = '%s'  " % chrom
    qry += "and start_position = %s  " % start
    rows = search_db(cursor, qry)
    store (cursor, header_fields, new_row)

    

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA", "FPPP", "GBM", "HNSC", "KICH" ,"KIRC","KIRP",
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    #check_strand_annotation(cursor, db_names)

    # how many cases do we have affected?
    total_cases = 0
    for db_name in db_names:
        print " ** ", db_name
        switch_to_db (cursor, db_name)

        # header fields:
        header_fields = []
        qry = "select column_name from information_schema.columns where "
        qry += "table_schema='%s' "  % db_name
        qry += "and table_name = 'somatic_mutations'"
        rows = search_db(cursor, qry)
        for  row in rows:
            header_fields.append(row[0])

 
        qry  = "select distinct  tumor_sample_barcode from somatic_mutations "
        rows = search_db(cursor, qry)
        per_db_cases = 0
        if not rows: continue
        for  row in rows:
            tbarcode = row[0]
            # this is wrong - not taking different cromosomes into accoutn
            qry =  'select distinct chromosome, start_position, count(*) from somatic_mutations '
            qry += "where tumor_sample_barcode = '%s' " % tbarcode
            qry += "group by start_position having count(*) > 1 "
            rows2 = search_db (cursor, qry)
            if not rows2: continue

            for row2 in rows2:
                [chrom, start, count] = row2
                print chrom, start, count
                #exit(1)
                #cleanup (cursor, header_fields, tbarcode, chrom, start)
                per_db_cases += 1


        print "cases: ", per_db_cases
        print
    
        total_cases += per_db_cases

    print
    print "total cases: ", total_cases
    print

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
