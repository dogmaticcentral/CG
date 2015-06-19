#!/usr/bin/python -u

# to speed up things make index:  (int 01_maf_somatic_mutations_table.sql)
# create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position)
# note: I am assuming then that no tumot will have the exact same mutation in both alleles, for the simple reason that I do nto see
# how would the thwo entried in the database then be distinguished from one another
# (rather if tumor_seq_allele1 == tumor_seq_allele2 != match_norm_seq_allele1  and tumor_seq_allele2 != match_norm_seq_allele2
# then I have somehting like that)


import os
from   tcga_utils.mysql   import  *
import commands
from os import listdir
from time import time


#########################################
def check_barcode(cursor, barcode):

    # are there any other annotations (usually indicating something is off with the sample
    switch_to_db(cursor, 'tcga_meta')
    qry = "select * from annotation where item_barcode='%s'" % barcode
    rows = search_db(cursor, qry)
    if rows:
        print "annotation found:"
        print rows
        return False

    return True

#########################################
def store (cursor, header_fields, fields):

    if len(fields) != len(header_fields): return

    all_fields = {}
    for i in range( len(header_fields)):
        field = fields[i]
        header = header_fields[i]
        all_fields[header] = field

    ok = store_without_checking (cursor, 'rnaseq_rpkm', all_fields)

    if not ok:
        print 'store failure:'
        print all_fields
        exit(1)

#########################################
def barcode2short(barcode):
    # barcode short            # the elemnts of the barcode are
    # project - tissue source site (TSS)  - participant -
    # source.vial - portion.analyte  - plate - (sequencing or characterization center)
    elements = barcode.split('-')
    source  = elements[3][:-1]
    # source can signa additional or metastatic tumors from the same patient
    # to keep our life simple we'll just stick to primary tumors
    # indicated by source code 01, 03, 08, or 09
    sample_barcode_short = '-'.join(elements[1:3] + [source]) # get rid of the 'vial' character

    return [sample_barcode_short, int(source)]

#########################################
def load_expression_file (cursor, db_name, header_fields,  barcode_in, expr_file):

    switch_to_db (cursor, db_name)
    fields = expr_file.split('/')[-1].split('.')
    barcode = fields[1]
    if barcode_in != barcode:
        print "barcode-file mismatch:"
        print barcode
        print expr_file
        exit(1)
    # this happens to work for files deposited by UNC, and there are the only ones that I have - not sure
    # if this is going to work in  general
    experiment_id = fields[0]
    [sample_barcode_short, source_code] = barcode2short(barcode)

    inff = open(expr_file, "r")
    start = time ()
    for line in inff.readlines( )[1:]:
        if line[0] == '?': continue
        line = line.rstrip()
        fields = line.split ('\t')
        if len(fields) != 4: continue # I don't know what this is in that case
        fields_clean = [x.replace("'", '') for x in fields]
        symbol = fields_clean[0].split('|')[0]
        rpkm = float (fields_clean [-1])

        store (cursor, header_fields, [symbol, barcode, sample_barcode_short, rpkm, source_code, experiment_id])

    inff.close()
    print barcode, "done in %5.2fs" % (time() - start)

#########################################
def check_duplicates_and_flagged_barcodes (cursor, parent_dirs, files_in_data_dir):
    data_files_per_barcode = {}
    for parent_dir in parent_dirs:
        # the magical formula in check_barcode gives me the barcode from the name of the file
        filtered_files = [x for x in files_in_data_dir[parent_dir] if check_barcode(cursor, x.split('.')[1]) ]
        print "parent dir:", parent_dir, "number of files:", len(filtered_files)
        for data_file in filtered_files:
            fields  = data_file.split('.')
            barcode = fields[1]
            if not data_files_per_barcode.has_key(barcode):
                data_files_per_barcode[barcode] = [parent_dir + "/" + data_file]
            else:
                data_files_per_barcode[barcode].append(parent_dir + "/" + data_file)
    return data_files_per_barcode

#########################################
def average_and_load (cursor, db_name, header_fields, barcode, files):

    switch_to_db (cursor, db_name)
    [sample_barcode_short, source_code] = barcode2short(barcode)
    start = time ()
    values = {}
    experiment_id = ""
    for expr_file in files:
        fields = expr_file.split('/')[-1].split('.')
        if experiment_id: experiment_id += ","
        experiment_id += fields[0]
        inff = open(expr_file, "r")

        for line in inff.readlines()[1:]:
            if line[0] == '?': continue
            line = line.rstrip()
            fields = line.split ('\t')
            if len(fields) != 4: continue # I don't know what this is in that case
            fields_clean = [x.replace("'", '') for x in fields]
            # TODO add symbol name resolution here
            symbol = fields_clean[0].split('|')[0]
            rpkm = float (fields_clean [-1])
            if not values.has_key(symbol): values[symbol] = []
            values[symbol].append(rpkm)
        inff.close()

    for symbol, vals in values.iteritems():
        avg = 0.0
        if len(vals):
            for v in vals: avg += v
            avg /= len(vals)
        #print "%10s  %5.2f" % (symbol, avg), vals
        store (cursor, header_fields, [symbol, barcode, sample_barcode_short, avg, source_code, experiment_id])
    print barcode, "(with averaging) done in %5.2fs" % (time() - start)



#########################################
def process_data_set (cursor, db_name, data_files_per_barcode):

    header_fields = ['symbol', 'tumor_sample_barcode', 'sample_barcode_short', 'rpkm', 'source_code', 'experiment_id']

    for barcode, files in data_files_per_barcode.iteritems():
        # the magical formula in check_barcode gives me the barcode from the name of the file
        if len(files) == 1:
            load_expression_file (cursor, db_name, header_fields, barcode, files[0])
        else:
            # this is much cheaper (mysql-wise) than storing both
            average_and_load (cursor, db_name, header_fields,  barcode, files)

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","REA","UCEC"]

    print "NOTE <<<<>>>>> duplicate checking is disabled on the db level - this"
    print "has to be done on an empyt database !!!!!"
    exit(1)

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

        parent_dirs = []
        files_in_data_dir = {}
        # the exact same files (to the number) appear here not under different revision,
        # but under different set number, ie.
        #_BRCA.IlluminaHiSeq_RNASeq.Level_3.1.2.0 and _BRCA.IlluminaHiSeq_RNASeq.Level_3.7.0.0
        # both contain UNCID_602432.TCGA-A8-A091-01A-11R-A00Z-07.111011_UNC12-SN629_0149_AB062VABXX.8.trimmed.annotated.gene.quantification.txt
        # I have to make some decision and move on - I decide to keep the file when it appears for the first time
        seen = {}
        for parent_dir, subdir_list, file_list in os.walk(db_dir):
            data_files = filter (lambda x: x[-3:] == 'txt', file_list)
            df_clean = []
            for df in data_files:
               if seen.has_key(df):
                   pass
               else:
                   df_clean.append(df)
                   seen[df] = True
            if df_clean:
                parent_dirs.append(parent_dir)
                files_in_data_dir[parent_dir] = df_clean

        if not parent_dirs:
            print "no data found for db_name (?)"
        else:
            data_files_per_barcode = check_duplicates_and_flagged_barcodes (cursor, parent_dirs, files_in_data_dir)
            process_data_set (cursor, db_name, data_files_per_barcode)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
