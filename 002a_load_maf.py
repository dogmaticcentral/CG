#!/usr/bin/python

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
def store (cursor, header_fields, fields, usable_field):
    
    fixed_fields  = {}
    update_fields = {}

    # for some fiels, the length of the header row and the 
    # data rows is not equal
    # there isn't much I can do about it:
    shorter =  len(header_fields) 
    if  len(fields) <  shorter: 
        shorter =  len(fields) 
    
    for i in range(shorter ):
        if not usable_field[i]: continue
        field = fields[i]
        header = header_fields[i]
        # if we first bump into data set without aa_code info, we just store it;
        # in the following script we'll resolve it and keep the more  more informative
        # field from each entry
        if (header in ['tumor_sample_barcode', 'chromosome', 'aa_code', 'start_position'] ):
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
def load_maf (cursor, db_name, required_fields, maffile):
    inff = open(maffile, "r")
    line_ct = 0
    missing_fields = []
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
                elif header_fields[i] in ['cdna_change', 'chromchange', 
                                          'c_position_wu', 'c_position']:
                    header_fields[i] = 'cdna_change'

            missing_fields =  list (Set(required_fields) - Set(header_fields) )
            # difference: fields in usable_fields, but not in header_fields
            if  len(missing_fields) > 0:
                print "   missing fields: ", Set(required_fields) - Set(header_fields) 
            else:
                print "   no fields missing"
            header_fields += missing_fields
            # note: usable_field is a list of booleans
            usable_field  = map (lambda x: x in required_fields, header_fields)

        else:
            fields_clean = [x.replace("'", '') for x in fields]
            for miss in missing_fields:
                fields_clean.append("missing")
            store (cursor, header_fields, fields_clean, usable_field)
    inff.close()



#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA", "FPPP", "GBM", "HNSC", "KICH" ,"KIRC","KIRP",
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    for db_name in db_names:
        # check db exists
        #qry = "show databases like '%s'" % db_name
        #rows = search_db(cursor, qry)
        #if not rows:
        #    print db_name, "not found"
        #    exit(1)
        

        print " ** ", db_name
        switch_to_db (cursor, db_name)

        table = 'somatic_mutations'

        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name
            exit(1)

        # if there is a nonexistent field from the ones that we require, drop the whole dataset
        # (we want only the cases where the mutation is traceable back to the nucleotide position in cDNA,
        # and for several 'blacklisted' datasets this is not the case)
        required_fields = find_existing_fields(cursor, db_name)

        db_dir  = '/Users/ivana/databases/TCGA/'+db_name
        ret       = commands.getoutput('find ' + db_dir + ' -name "*.maf"')
        maf_files = ret.split('\n')
        for maffile in maf_files:
            print '\t loading:', maffile
            load_maf (cursor, db_name, required_fields, maffile)

 
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
