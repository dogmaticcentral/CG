#!/usr/bin/python

import os.path
import commands
from tcga_utils.mysql  import  *
from tcga_utils.utils  import  get_expected_fields, process_header_line
from time import time
from random import random


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

    rows   = search_db (cursor, qry)
    # if there is a return,  something went wrong
    if (rows):
        search_db (cursor, qry, verbose=True)
        exit(1) # exit, bcs we should not be here

#########################################
def one_allele_normal(field) :
    tumor  = {field['tumor_seq_allele1'], field['tumor_seq_allele2']} # this is a set
    normal =  {field['match_norm_seq_allele1'], field['match_norm_seq_allele2']}
    return len(tumor&normal)>0

#########################################
def is_useful(fields, header):

    non_info  = ['missing', '', '.', '-']
    return fields != None and fields.has_key(header) and fields[header] != None and  not fields[header].replace(" ", "") in non_info

#########################################
def conflict_annotation_updated (existing_fields, conflict_comment, new_id):
    conflict_annotation = existing_fields['conflict']
    if not conflict_annotation:
        conflict_annotation = ""
    elif len(conflict_annotation) > 0:
        conflict_annotation += "; "
    conflict_annotation += "%s with %d" % (conflict_comment, new_id)
    return conflict_annotation

#########################################
def update_conflict_field (cursor, table, existing_row_id, existing_fields, new_id, conflict_comment):
    conflict_annotation = conflict_annotation_updated (existing_fields, conflict_comment, new_id)
    update_fields = {'conflict': conflict_annotation}
    existing_fields['conflict'] = conflict_annotation
    update_db (cursor, table, existing_row_id, update_fields)

#########################################
def resolve_duplicate (cursor, table, expected_fields, existing_rows, new_fields):

    this_is_a_duplicate = False
    diagnostics = {}

    existing_fields_by_database_id  = dict( zip (map (lambda x: x[0], existing_rows),  map (lambda x: make_named_fields(expected_fields, x[1:]), existing_rows) ))

    # try to diagnose how come we have multiple reports for a mutation starting at a given position
    for db_id, existing_fields in  existing_fields_by_database_id.iteritems():
        # the end position the same?
        diagnostics[db_id] = ""
        if existing_fields['end_position'] == new_fields['end_position']:
            # are all alleles the same?
            all_alleles_the_same = True
            empty_alleles = {'existing_fields': 0, 'new_fields': 0}
            for allele  in ['reference_allele', 'tumor_seq_allele1', 'tumor_seq_allele2',
                                'match_norm_seq_allele1', 'match_norm_seq_allele2']:
                all_alleles_the_same = all_alleles_the_same and existing_fields[allele] == new_fields[allele]
                # are any of the entries empty?
                if existing_fields[allele]==None or len(existing_fields[allele])==0 :
                    empty_alleles['existing_fields'] += 1
                if new_fields[allele]==None or len(new_fields[allele])==0 :
                    empty_alleles['new_fields'] += 1

            if all_alleles_the_same:
                if existing_fields['variant_classification'] != new_fields['variant_classification']:
                    diagnostics[db_id] = "conflict: different variant classification"
                else:
                    # do both entries have the same aa_change info
                    if new_fields.has_key('aa_change'):
                        if existing_fields['aa_change'] == new_fields['aa_change']:
                            #this is an exact duplicate, do nothing
                            diagnostics[db_id] = "move on"
                            this_is_a_duplicate = True
                            break
                        elif is_useful (existing_fields, 'aa_change') and not is_useful (new_fields, 'aa_change'):
                            diagnostics[db_id] = "move on: the old has the aa info"
                            this_is_a_duplicate = True
                            break
                        elif not is_useful (existing_fields, 'aa_change') and is_useful (new_fields, 'aa_change'):
                            diagnostics[db_id] = "use new: it has the aa info"
                        else:
                            diagnostics[db_id] = "conflict: different aa"
                            # should I be resolving it at this place?
                    else:
                        diagnostics[db_id] = "conflict: of unclear origin"

            else: # some alleles are different
                diagnostics[db_id] = "different alleles;  empty existing: %d, empty new: %d      " % (empty_alleles['existing_fields'], empty_alleles['new_fields'])
                if empty_alleles['existing_fields'] > empty_alleles['new_fields']:
                    # the new entry has more info
                    diagnostics[db_id] += "use new"
                else:
                    diagnostics[db_id] += "move on: keep old"
                    this_is_a_duplicate = True
                    break

        else: # the end position is not the same
            diagnostics[db_id] = "end positions different   *%s*   *%s* " % ( existing_fields['end_position'], new_fields['end_position'])
            # how could that happen?
            # one of the possibilities (gosh, will I have to go through all of them?)
            # is that a frameshift mutation is interpreted differently
            # hard to find a robust solution
            if existing_fields['variant_classification'] == 'frame_shift_del' and  one_allele_normal(existing_fields) and one_allele_normal(new_fields):
                    # do nothing, this is a different interpretation of the same mutation
                    diagnostics[db_id] = "move on"
                    this_is_a_duplicate = True
                    break

            elif not one_allele_normal(existing_fields) and  not one_allele_normal(new_fields):
                    # store, this is possibly compound  heterozygous
                    diagnostics[db_id] = "compound heterozygous"
            else:
                # I don't know what this is
                diagnostics[db_id] = "conflict: different length"


    if False and new_fields['start_position'] == 35043650:
        print
        print "="*20
        print
        for header in new_fields.keys():
            print "%30s    [ %s ] " % (header, new_fields[header]),
            for existing_fields in existing_fields_by_database_id.values():
                print " [ %s ] " % (existing_fields[header]),
            print
        print
        print "this is a duplicate ", this_is_a_duplicate
        print
        if this_is_a_duplicate: return "duplicate"
        for db_id in existing_fields_by_database_id.keys():
            print db_id, "    ", diagnostics[db_id]
        print

    # if it is a duplicate of any of the existing entries, we do not have to worry about it any more
    if this_is_a_duplicate: return "duplicate"

    for existing_row_id, existing_fields in existing_fields_by_database_id.iteritems():
        if 'use new' in diagnostics[existing_row_id]:
            update_fields = new_fields
            update_db (cursor, table, existing_row_id, update_fields)
            return "to replace"  # better info than what we already have
        if 'compound' in diagnostics[existing_row_id]:
            new_id = insert_into_db (cursor, table, new_fields)
            update_conflict_field (cursor, table, existing_row_id, existing_fields, new_id, "compound")
            update_conflict_field (cursor, table, new_id, new_fields, existing_row_id, "compound")
            return "compound"  # store both without alarm - I'm not sure that this ever happens

    # if we did not return until this point, we have a conflict - the question is with how many rows
    new_id = -1
    for existing_row_id, existing_fields in existing_fields_by_database_id.iteritems():
        if 'conflict' in diagnostics[db_id]:
            if new_id<0 :   new_id = insert_into_db (cursor, table, new_fields)
            update_conflict_field (cursor, table, existing_row_id, existing_fields, new_id, "unresolved")
            # the "new" raw already also became existing by this point
            update_conflict_field (cursor, table, new_id, new_fields, existing_row_id, "unresolved")

    return "conflict"

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
def make_named_fields (header_fields, fields, expected_fields = None):

    named_fields = {}

    if len(header_fields) != len(fields) :
        print "##################################"
        print "fields length mismatch (?)" # it should have been solved by this point
        print len(header_fields), len(fields)
        exit(1) # header field mismatch

    for i in range( len(header_fields)):
        header = header_fields[i]
        if expected_fields and not header in expected_fields: continue
        field = fields[i]
        named_fields[header] = field

    return named_fields

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


    if not existing_rows: #this is the first time we a mutation in this place in this patient
        insert_into_db(cursor, table, new_fields)
        return "new"
    else:
        # do something about possible duplicates
        return resolve_duplicate (cursor, table, expected_fields, existing_rows, new_fields)

    return ""

#########################################
def field_cleanup(maf_header_fields, sample_barcode_short, meta_id, maf_fields):

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

    # also in addition to the fields in the original maf file, we are adding a pointer to the maffile itself
    clean_fields.append(meta_id)
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
            clean_fields[cdna_change_field_idx] = cdna_change_field_entry.replace ("c.", "")
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
def load_maf (cursor, db_name, table, maffile, meta_id, stats):

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

    expected_fields = get_expected_fields(cursor, db_name, table)
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

        clean_fields = field_cleanup(maf_header_fields, sample_barcode_short, meta_id, maf_fields)
        retval = store(cursor, table, expected_fields,
              maf_header_fields + ['sample_barcode_short', 'meta_info_index', 'conflict'], clean_fields)
        if not retval in stats.keys(): stats[retval] = 0
        stats[retval] += 1

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
                 "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",  "MESO", "OV",   "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    #db_names  = ["ACC"]

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

        qry = "select * from mutations_meta"
        rows = search_db(cursor, qry)
        if not rows:
            print "no meta info found"
            continue

        db_dir  = '/mnt/databases/TCGA'
        if not  os.path.isdir(db_dir):
            print "directory " + db_dir + " not found"
            exit(1) # TCGA db dir not found

        maf_file = {}
        for row in rows:
            [meta_id, file_name, quality_check, assembly, diagnostics] = row
            if quality_check=="fail": continue
            maf_file[meta_id] = "/".join([db_dir, db_name, "Somatic_Mutations", file_name])

        stats = {}
        for meta_id, maffile in maf_file.iteritems():
            load_maf (cursor, db_name, table, maffile, meta_id, stats)
        print
        for type, count in stats.iteritems():
            print type, count

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
