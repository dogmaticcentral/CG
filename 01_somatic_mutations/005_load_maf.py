#!/usr/bin/python

#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2016 Ivana Mihalek.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#

import os.path
import re, commands
from tcga_utils.utils  import  *
from time import time

dorky = re.compile('(\-*\d+)([ACGT]+)>([ACGT]+)')
#########################################

def clean_cdna_change_annotation(old_annot):
    new_annot = old_annot.replace ("c.", "")
    if '>' in new_annot:
        match_return = re.match(dorky,new_annot)
        if not match_return:
            # I have no idea what this is
            # at least we got rid of the 'c.' crap
            pass
        else:
            new_annot = "%s%s%s" %( match_return.group(2),  match_return.group(1),  match_return.group(3))
    return new_annot

#########################################
def update_db (cursor, table, row_id, update_fields):

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
    rows = search_db (cursor, qry)

    # if there is a return,  something went wrong
    if (rows):
        search_db (cursor, qry, verbose=True)
        exit(1) # exit, bcs we should not be here


#########################################
def insert_into_db (cursor, table, fields):

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

    qry = "select last_insert_id()"
    rows = search_db (cursor, qry)
    if not rows:
        print "last insert id failure (?)"
        exit(1) # last insert id failure

    return int(rows[0][0])

#########################################
def store (cursor, table, expected_fields, maf_header_fields, new_row):

    sample_barcode_short = new_row[ maf_header_fields.index('sample_barcode_short')]
    chromosome           = new_row[ maf_header_fields.index('chromosome')]
    start_position       = new_row[ maf_header_fields.index('start_position')]

    qry  = "select * from %s  " % table
    qry += "where sample_barcode_short = '%s' " % sample_barcode_short
    qry += "and chromosome = '%s'  " % chromosome
    qry += "and start_position = %s  " % start_position
    existing_rows = search_db(cursor, qry)

    new_fields = make_named_fields(maf_header_fields, new_row, expected_fields)

    if not existing_rows: #this is the first time we see a mutation in this place in this patient
        insert_into_db(cursor, table, new_fields)
        return "new"
    else:
        # do something about possible duplicates
        print "duplicate!"
        #exit(1)
        #  apparently, they have cleaned up; this happens so rarely that I wil just ignore it
        # I have seen only one 2 cases in UCEC (and nowhere else)
        return "duplicate"
        #return resolve_duplicate (cursor, table, expected_fields, existing_rows, new_fields)


    return ""

#########################################
def field_cleanup(maf_header_fields, sample_barcode_short, maf_fields):

    number_of_header_fields = len(maf_header_fields)

    # TCGA is exhausting all possibilities here (row longer than the header):
    clean_fields = maf_fields[:number_of_header_fields]  # I do not know what you guys are anyway, so off you go
    clean_fields = [x.replace("'", '') for x in clean_fields]  # get rid of the quote marks
    clean_fields = [x.replace(" ", '') for x in clean_fields]  # get rid of the quote marks

    # is the number of the fields smaller than the number we are expecting from the header?
    # (yes, that can happen, anybody can submit anything in whichever fucking format they like in here)
    for i in range(len(clean_fields), number_of_header_fields):
        clean_fields.append("missing")

    clean_fields.append(sample_barcode_short)

    # there is an optional 'conflict field' which might get filled if there is one:
    clean_fields.append(None)
    # special: I want to be able to index on aa_change
    # so I limited the length of the aa_change field to 100 characters (if it does not know the
    # the length of the field, mysql refuses to index)
    # but sometimes people put large swath of sequence here; instead of chopping, replace with the mutation type
    if 'aa_change' in maf_header_fields:
        aa_change_field_idx = maf_header_fields.index('aa_change')
        aa_change_field_entry = clean_fields[aa_change_field_idx]
        if len(aa_change_field_entry) > 100:
            clean_fields[aa_change_field_idx] = clean_fields[maf_header_fields.index('variant_classification')]
        else:
            clean_fields[aa_change_field_idx] = aa_change_field_entry.replace ("p.", "")
    # the same for cdna
    if 'cdna_change' in maf_header_fields:
        cdna_change_field_idx = maf_header_fields.index('cdna_change')
        cdna_change_field_entry = clean_fields[cdna_change_field_idx]
        if len(cdna_change_field_entry) > 100:
            clean_fields[cdna_change_field_idx] = clean_fields[maf_header_fields.index('variant_classification')]
        else:
            clean_fields[cdna_change_field_idx] = clean_cdna_change_annotation(cdna_change_field_entry)
    # change to lowercase wherever possible
    for header in ['variant_classification', 'variant_type', 'verification_status', 'validation_status' , 'mutation_status']:
        index =  maf_header_fields.index(header)
        clean_fields[index] = clean_fields[index].lower()


    for header in ['start_position', 'end_position']:
        index =  maf_header_fields.index(header)
        clean_fields[index] = int(clean_fields[index])

    # one more thing, I hope it is the last
    chromosome_field = maf_header_fields.index('chromosome')
    if clean_fields[chromosome_field].upper() == "MT":
        clean_fields[chromosome_field] = "M"

    #tumor1_idx = maf_header_fields.index('tumor_seq_allele1')
    #tumor2_idx = maf_header_fields.index('tumor_seq_allele2')
    norm1_idx = maf_header_fields.index('match_norm_seq_allele1')
    norm2_idx = maf_header_fields.index('match_norm_seq_allele2')
    ref_idx   = maf_header_fields.index('reference_allele')


    var_type_idx  = maf_header_fields.index('variant_type')
    var_class_idx = maf_header_fields.index('variant_classification')
    if "del" in clean_fields[var_class_idx]:
        clean_fields[var_type_idx] = "del"
    if "ins" in clean_fields[var_class_idx]:
        clean_fields[var_type_idx] = "ins"

    # this all needs to be redone if they ever start top put in decent estimate for allales
    if clean_fields[var_type_idx] != "ins":
        if clean_fields[norm2_idx] == "-":
            clean_fields[norm2_idx] = ""
        if clean_fields[norm1_idx] == "-":
            clean_fields[norm1_idx] = ""


    return clean_fields

#########################################
def check_tumor_type(tumor_type_ids, maf_header_fields, maf_fields):

    # here is where we construct the short version of the barcode that identifies the sample
    # I am adding this one so I do not have to search the database by doing substring comparison
    tbarcode = maf_fields[maf_header_fields.index('tumor_sample_barcode')]
    # the elements of the barcode are
    # project - tissue source site (TSS)  - participant -
    # source.vial - portion.analyte  - plate - (sequencing or characterization center)
    elements = tbarcode.split('-')
    tumor_type_code = elements[3][:-1]
    # tumor_type_code can signal additional or metastatic tumors from the same patient
    # to keep our life simple we'll just stick to primary tumors
    # indicated by source code 01, 03, 08, or 09
    if not tumor_type_code in tumor_type_ids: return None

    sample_barcode_short = '-'.join(elements[1:3] + [tumor_type_code])  # get rid of the 'vial' character
    return sample_barcode_short

#########################################
def load_maf (cursor, db_name, table, maffile):

    if not os.path.isfile(maffile):
        print "not found: "
        print maffile
        exit(1)  # maffile not found
    cmd = "wc  -l  " + maffile
    nol = int(commands.getstatusoutput(cmd)[1].split()[0]) -1
    print "processing %d lines from %s " % (nol, maffile)

    tumor_type_ids = []
    if table == "somatic_mutations":
        tumor_type_ids = ['01', '03', '08', '09']
    elif table == "metastatic_mutations":
        tumor_type_ids = ['06']
    else:
        print "I don't know how to handle ", table, " sample type"
        exit(1) # unrecognized table name

    expected_fields   = get_expected_fields(cursor, db_name, table)
    maf_header_fields = process_header_line(maffile)

    inff = open(maffile, "r")
    line_ct = 0
    start = time()
    first = True
    tot_entries = 0
    type_queried = 0
    for line in inff:
        line_ct += 1
        if not line_ct%1000:
            print "\t processed %5d (%4.1f%%)  %8.2fs" % (line_ct, float(line_ct)/nol*100,  time()-start)
        if line.isspace(): continue
        if line[0]=='#': continue
        if first: # this is the header line
            first = False
            continue
        tot_entries += 1
        line = line.rstrip()
        maf_fields = line.split ('\t')

        sample_barcode_short = check_tumor_type(tumor_type_ids, maf_header_fields, maf_fields)
        if not sample_barcode_short: continue # this can happen if the tumor type is not the one we are looking for
        type_queried += 1

        clean_fields = field_cleanup(maf_header_fields, sample_barcode_short, maf_fields)
        retval = store(cursor, table, expected_fields,
              maf_header_fields + ['sample_barcode_short', 'conflict'], clean_fields)

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
        exit(1) # unknown sample type

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL",  "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",  "MESO", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
    db_names  = [ "UCEC", "UCS", "UVM"]

    #db_names = ["ACC"]

    for db_name in db_names:
        # check db exists
        qry = "show databases like '%s'" % db_name
        rows = search_db(cursor, qry)
        if not rows:
            print db_name, "not found"
            exit(1) # db not found

        print " ** ", db_name
        switch_to_db (cursor, db_name)
        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name

        db_dir = "/mnt/databases/TCGA/" + db_name + "/Somatic_Mutations"
        if not os.path.isdir(db_dir):
            print "directory " + db_dir + " not found"
            exit(1)

        mafs = filter(lambda x: x.endswith(".maf") and not 'mitochondria' in x, os.listdir(db_dir))
        maf_files =   ["/".join([db_dir,  maf]) for maf in mafs]
        for maffile in maf_files:
             load_maf (cursor, db_name, table, maffile)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
