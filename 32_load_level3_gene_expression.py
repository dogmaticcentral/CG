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
def store (cursor, header_fields, fields):
    
    fixed_fields  = {}
    update_fields = {}
    
    if (len(fields) != len(header_fields)): return
    for i in range( len(header_fields) ):
        field = fields[i]
        header = header_fields[i]
        if (header in ['symbol', 'tumor_sample_barcode', 'source'] ):
            fixed_fields[header] = field
        else:
            update_fields[header] = field
        
    ok = store_or_update (cursor, 'gene_expression', fixed_fields, update_fields)
    if not ok:
        print 'store failure:'
        print fixed_fields
        print update_fields
        exit(1)

 
#########################################
def load_expression_file (cursor, db_name,  expr_file, source):
    print expr_file
    inff = open(expr_file, "r")
    line_ct = 0
    tcga = ''
    header_fields = ['symbol', 'fold_change', 'tumor_sample_barcode', 'source']
    for line in inff:
        line = line.rstrip()
        fields = line.split ('\t')
        line_ct += 1
        if line_ct == 1:
            if not 'TCGA'  in fields[-1]:
                return
            print fields[-1] # this should be the barcode
            tcga = fields[-1]
        elif line_ct == 2:
            continue
        else:
            fields_clean = [x.replace("'", '') for x in fields]
            if len(fields_clean) != 2: continue # I don't know what this is in that case
            if fields_clean[1] == 'null': continue
            if not tcga:
                print 'sample_id not defined'
                exit(1)
                    
            fields_clean.append( tcga)
            fields_clean.append( source)
            store (cursor, header_fields, fields_clean)
    inff.close()
    


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["COAD",  # after this we go alphabetically
                 "ACC", "BLCA", "BRCA", "CESC",  "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]
    source = 'mdanderson'
    for db_name in db_names:
        # check db exists
        qry = "show databases like '%s'" % db_name
        rows = search_db(cursor, qry)
        if not rows:
            print db_name, "not found"
            exit(1)

        print " ** ", db_name
        switch_to_db (cursor, db_name)

        table = 'gene_expression'

        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name
            exit(1)

 
        db_dir  = '/Users/ivana/databases/TCGA/%s/Expression_Genes' % db_name
        cmd = 'ls ' + db_dir + "/%s*/*.txt" % source
        ret = commands.getoutput(cmd)
        if 'No such' in ret: 
            print ret
            continue
        if not ret: continue
        expression_files = ret.split('\n')
        for expr_file in expression_files:
            print '\t loading:', expr_file
            load_expression_file (cursor, db_name,  expr_file, source)
        exit(1)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
