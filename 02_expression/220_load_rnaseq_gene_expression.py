#!/usr/bin/python -u

# to speed up things make index:  (int 01_maf_somatic_mutations_table.sql)
# create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position)
# note: I am assuming then that no tumot will have the exact same mutation in both alleles, for the simple reason that I do nto see
# how would the thwo entried in the database then be distinguished from one another
# (rather if tumor_seq_allele1 == tumor_seq_allele2 != match_norm_seq_allele1  and tumor_seq_allele2 != match_norm_seq_allele2
# then I have somehting like that)


import MySQLdb
from   tcga_utils.mysql   import  *
import commands
from os import listdir
from time import time

#########################################
def store (cursor, header_fields, fields):
    
    fixed_fields  = {}
    update_fields = {}
    
    if (len(fields) != len(header_fields)): return
    for i in range( len(header_fields) ):
        field = fields[i]
        header = header_fields[i]
        if (header in ['symbol', 'sample_id', 'experiment_id'] ):
            fixed_fields[header] = field
        else:
            update_fields[header] = field
        
    ok = store_or_update (cursor, 'rnaseq_rpkm', fixed_fields, update_fields)
    if not ok:
        print 'store failure:'
        print fixed_fields
        print update_fields
        exit(1)

 
#########################################
def load_expression_file (cursor, db_name,  expr_file):
    #print expr_file
    inff = open(expr_file, "r")
    fields = expr_file.split('/')[-1].split('.')
    barcode = fields[1]
    # this happens to work for files deposited by UNC, and there are the only ones that I have - not sure
    # if this is going to work in  general
    experiment_id = fields[0]

    fields =  barcode.split ('-')
    source_code = int(fields[3][:2])
    patient = '-'.join(fields[:3])
    # paired sample?
    if (source_code== 11): # this is normal, it may correspond to several tumor samples
        paired_barcode = ""
    else:
        directory = '/'.join( expr_file.split('/')[:-1] )
        normal_barcodes = []
        for file in listdir(directory):
            barcode2 = file.split('.')[1]
            if patient + "-11" in barcode2:
                normal_barcodes.append(barcode2)
        # there should be one and only one normal barcode
        if len(normal_barcodes) == 0:
            print "no normal match found for ", barcode
            exit(1)
        elif len(normal_barcodes) > 1:
            print "multiple normal found for ", barcode, "- ", normal_barcodes
            exit(1)
        else:
            paired_barcode = normal_barcodes[0]


    header_fields = ['symbol', 'sample_id', 'paired_sample_id', 'rpkm', 'source_code', 'experiment_id']
    for line in inff.readlines( )[1:]:
        if line[0] == '?': continue
        line = line.rstrip()
        fields = line.split ('\t')
        if len(fields) != 4: continue # I don't know what this is in that case
        fields_clean = [x.replace("'", '') for x in fields]
        symbol = fields_clean[0].split('|')[0]
        rpkm = float (fields_clean [-1])

        store (cursor, header_fields, [symbol, barcode, paired_barcode, rpkm, source_code, experiment_id])
    inff.close()
    


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","UCEC"]


    for db_name in db_names:
        # check db exists
        qry = "show databases like '%s'" % db_name
        rows = search_db(cursor, qry)
        if not rows:
            print db_name, "not found"
            exit(1)

        print " ** ", db_name
        switch_to_db (cursor, db_name)

        table = 'rnaseq_rpkm'

        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name
            exit(1)

 
        db_dir  = '/Users/ivana/databases/TCGA/%s/Expression_Genes' % db_name
        cmd = 'ls ' + db_dir + "/*/*.txt"
        ret = commands.getoutput(cmd)
        if 'No such' in ret: 
            print ret
            continue
        if not ret: continue
        expression_files = ret.split('\n')
        for expr_file in expression_files:
            print '\t loading:', expr_file.split('/')[-1].split('.')[1],
            #cmd = 'head ' + expr_file
            #ret = commands.getoutput(cmd)
            #print ret
            #print
            start = time()
            load_expression_file (cursor, db_name,  expr_file)
            print "\t\t done in %.2f secs" % (time() - start)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
