#!/usr/bin/python


# damage control: 
# Originally I had the stupid idea (or, I'd like to think that I just wasn't thinking)
# to do the test for existing entries like this:
#        if (header in ['tumor_sample_barcode', 'chromosome', 'strand', 'start_position'] ):
#            fixed_fields[header] = field
# However, this would give me duplicates, because some submitters may call the strand '+', and some may call it '+1'
# The complete list of symbols used to denote strand:  ['+', '+1', '-1', '1', '-', ''].
# But then some put in empty string, which might be right: if a mutation exists, it will be present in both
# strands - should not have used strand as the criterion for uniqueness at all. I should have just dropped it.

# But. But then some good came out of it, bcs as it turns out, some centers deposit the info about the involved AA mutation,
# which makes my life significantly easier. So the idea now is the following. Do the database filling in two passes.
# In the first pass just stuff everything in (this was done in 002a_load_maf.py). Except, in distinction to my
# original stupid approach, now I am determining the uniquenes by the following fields:
# 'sample_barcode_short', 'chromosome',  'start_position', 'aa_change'
# Then, here, I check the duplicates, and
# combine the provided info.

# There is a serious issue that I am shoving under the carpet here:
# what if the two sets of somatic mutations are significantly different?
# It shouldn't be my call which one is different. Therefore I should probably throw away the data.
# I am solving it by pretending that the problem does not exist: the wo sets complement each other, and the difference
# is negligible.


import MySQLdb
from sets import Set
from   tcga_utils.mysql   import  *
import commands


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    table_name = 'somatic_mutations'
    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL",  "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",  "MESO", "OV",   "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    # how many cases do we have affected?
    total_cases = 0
    for db_name in db_names:
        print " ** ", db_name
        switch_to_db (cursor, db_name)

        # header fields:
        header_fields = []
        qry = "select column_name from information_schema.columns where "
        qry += "table_schema='%s' "  % db_name
        qry += "and table_name = '%s'" % table_name
        rows = search_db(cursor, qry)
        for  row in rows:
            header_fields.append(row[0])

 
        qry  = "select distinct  sample_barcode_short from %s " % table_name
        rows = search_db(cursor, qry)
        per_db_cases = 0
        if not rows: continue
        for  row in rows:
            short_barcode = row[0]
            # this is wrong - not taking different cromosomes into accoutn
            qry =  'select chromosome, start_position, count(*) from %s ' % table_name
            qry += "where sample_barcode_short = '%s' " % short_barcode
            qry += "group by chromosome, start_position having count(*) > 1 "
            rows2 = search_db (cursor, qry)
            if not rows2: continue

            for row2 in rows2:
                [chrom, start, count] = row2
                print ">>>>>> ", short_barcode, chrom, start, count
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
