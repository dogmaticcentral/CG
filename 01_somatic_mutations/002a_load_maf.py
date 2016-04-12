#!/usr/bin/python

import os.path
from sets import Set
from tcga_utils.mysql   import  *


#########################################
def find_existing_fields(cursor, db_name, table):

    qry  = "describe %s " % table
    rows = search_db (cursor, qry)
    if not rows:
        print "%s not found in  %s" % (table, db_name)
        exit(1)
    # the first header field is 'id' - this is entry id in the table
    existing_fields = [row[0] for row in rows[1:]]
    return existing_fields

#########################################
def insert_into_db (cursor, fields, table):

    qry = "insert into %s "  % table
    qry += "("
    first = True
    for field in fields.keys(): # again will have to check for the type here
        if (not first):
            qry += ", "
        qry += field
        first = False
    qry += ")"

    qry += " values "
    qry += "("
    first = True
    for value in fields.values(): # again will have to check for the type here
        if (not first):
            qry += ", "
        if  value is None:
            qry += " null "
        elif type(value) is int:
            qry += " %d" % value
        elif type(value) is float:
            qry += " %f" % value
        else:
            qry += " \'%s\'" % value
        first = False
    qry += ")"

    rows   = search_db (cursor, qry)
    # if there is a return something went wrong
    if (rows):
        search_db (cursor, qry, verbose=True)
        exit(1) # exit, bcs we should not be here


#########################################
def update_db (cursor, row_id, update_fields, table):

    qry = "update %s set " % table
    first = True
    for field, value in update_fields.iteritems():
        if (not first):
            qry += ", "
        qry += " %s = " % field
        if  value is None:
            qry += " null "
        elif type(value) is int:
            qry += " %d" % value
        else:
            qry += " \'%s\'" % value
        first = False

    qry += " where id = %d" % int(row_id)

    rows   = search_db (cursor, qry)
    # if there is a return,  something went wrong
    if (rows):
        search_db (cursor, qry, verbose=True)
        exit(1) # exit, bcs we should not be here


#########################################
def make_named_fields (header_fields, fields, usable_field):

    named_fields = {}

    # in some cases, the length of the header row and the
    # data rows is not equal
    # there isn't much I can do about it:
    shorter =  len(header_fields)
    if  len(fields) <  shorter:
        shorter =  len(fields)

    for i in range(shorter ):
        if not usable_field[i]: continue
        field = fields[i]
        header = header_fields[i]
        named_fields[header] = field

    return named_fields

#########################################
def resolve_duplicate (cursor, db_fields, db_row, new_fields, table):

    non_info = ['missing', '', '.', '-']
    row_id = db_row[0]
    db_row = db_row[1:]

    no_cols = len(db_fields)
    if len (new_fields) != no_cols:
        print "oink ?!"
        exit(1)

    update_fields = {}
    for i in range(no_cols):
        header = db_fields[i]
        field  = db_row[i]
        if  not field in non_info: continue
        new_field = new_fields[header]
        if  new_field  in non_info: continue
        update_fields[header] = new_field
        #print "replacing |  ", field, " |  with ",  new_field

    if update_fields:
        update_db (cursor, row_id, update_fields, table)
        #exit(1)
    else:
        #print "no update"
        pass
    return


#########################################
def store (cursor, required_fields, header_fields, fields, usable_field, table):


    sample_barcode_short = fields[ header_fields.index('sample_barcode_short')]
    chromosome           = fields[ header_fields.index('chromosome')]
    start_position       = fields[ header_fields.index('start_position')]

    qry  = "select * from %s  " % table
    qry += "where sample_barcode_short = '%s' " % sample_barcode_short
    qry += "and chromosome = '%s'  " % chromosome
    qry += "and start_position = %s  " % start_position
    rows = search_db(cursor, qry)

    #print "\t ", sample_barcode_short, chromosome, start_position,

    if not rows:
        # store
        #print "storing"
        insert_into_db(cursor, make_named_fields(header_fields, fields, usable_field), table)

    elif len(rows) > 1:
        print "More than two rows for ", sample_barcode_short, chromosome, start_position
        exit(1)

    else:
        # do something about duplicates
        # required_fields are the actual header in the db table
        resolve_duplicate (cursor, required_fields, rows[0], make_named_fields(header_fields, fields, usable_field), table)

    return


#########################################
import commands
from time import time

def load_maf (cursor, db_name, required_fields, maffile, table):

    if not os.path.isfile(maffile):
        print "not found: "
        print maffile 
        exit(1)
    cmd = "wc  -l  " + maffile
    nol = int(commands.getstatusoutput(cmd)[1].split()[0]) -1
    print "\t number of entries:  %d " % nol


    if table == "somatic_mutations":
        source_codes = ['01', '03', '08', '09']
    elif table == "metastatic_mutations":
        source_codes = ['06']
    else:
        print "I don't know how to handle ", table, " sample type"
        exit(1)

    inff = open(maffile, "r")
    missing_fields = []
    header_fields  = []
    expected_number_of_fields  = 0
    line_ct = 0
    start = time()
    first = True
    for line in inff:
        line_ct += 1
        if not line_ct%1000:
            print "\t processed %5d (%4.1f%%)  %8.2fs" % (line_ct, float(line_ct)/nol*100,  time()-start)
        if line.isspace(): continue
        if line[0]=='#': continue
        line = line.rstrip()
        fields = line.split ('\t')
        if first:
            first = False
            expected_number_of_fields = len(fields)
            header_fields = [x.lower() for x in fields]
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
            # is the number of the fields smaller than the number we are expecting from the header?
            # (yes, that can happen, anybody can submit anything in whichever fucking format they like in here)
            for i in range(len(fields), expected_number_of_fields):
                fields.append('')
            # TCGA is exhausting all possibilities here (row longer than the header):
            if len(fields) > expected_number_of_fields:
                fields = fields [:expected_number_of_fields] # I do not know what you guys are anyway, so off you go
            fields_clean = [x.replace("'", '') for x in fields]
            # here is where we construct the short version of the barcode that identifies the sample
            tbarcode = fields[ header_fields.index('tumor_sample_barcode')]
            # the elemnts of the barcode are
            # project - tissue source site (TSS)  - participant -
            # source.vial - portion.analyte  - plate - (sequencing or characterization center)
            elements = tbarcode.split('-')
            source  = elements[3][:-1]
            # source can signa additional or metastatic tumors from the same patient
            # to keep our life simple we'll just stick to primary tumors
            # indicated by source code 01, 03, 08, or 09
            if not source in source_codes: continue
            sample_barcode_short = '-'.join(elements[1:3] + [source]) # get rid of the 'vial' character


            fields_clean.append(sample_barcode_short)
            for miss in missing_fields:
                fields_clean.append("missing")
            # special: I want to be able to index on aa_change (just because of the magical combo above),
            # so I limited the length of the aa_change field to 100 characters (if it does not know the
            # the length of the field, mysql refuses to index)
            # but sometimes people put large swath of sequence here; instead of chopping, replace with the mutation type
            aa_change_field = header_fields.index('aa_change')
            if  len(fields_clean[aa_change_field]) > 100:
                fields_clean[aa_change_field] = fields_clean[  header_fields.index('variant_classification')]

            # one more thing, I hope it is the last
            chromosome_field = header_fields.index('chromosome')
            if fields_clean[chromosome_field].upper() == "MT":
                fields_clean[chromosome_field] = "M"

            store (cursor, required_fields, header_fields, fields_clean, usable_field, table)
    inff.close()



#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    sample_type = "primary"


    if sample_type == "primary":
        table = 'somatic_mutations'
    elif sample_type == "metastatic":
        table = 'metastatic_mutations'
    else:
        print "I don't know how to hadndle ", sample_type, " sample types"
        exit(1)


    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL",  "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",  "MESO", "OV",   "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]


    for db_name in db_names:
        # check db exists
        qry = "show databases like '%s'" % db_name
        rows = search_db(cursor, qry)
        if not rows:
            print db_name, "not found"
            exit(1)
        

        print " ** ", db_name
        switch_to_db (cursor, db_name)



        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name
            continue

        # if there is a nonexistent field from the ones that we require, drop the whole dataset
        # (we want only the cases where the mutation is traceable back to the nucleotide position in cDNA,
        # and for several 'blacklisted' datasets this is not the case)
        required_fields = find_existing_fields(cursor, db_name, table)

        db_dir  = '/mnt/databases/TCGA/' + db_name
        if not  os.path.isdir(db_dir):
            print "directory " + db_dir + " not found"
            exit(1)
        ret       = commands.getoutput('find ' + db_dir + ' -name "*.maf"')
        maf_files = ret.split('\n')
        for maffile in maf_files:
            print '\t loading:', maffile.split('/')[-1]
            load_maf (cursor, db_name, required_fields, maffile, table)

 
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
