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
from tcga_utils.mysql  import  *
from tcga_utils.utils  import  get_expected_fields, process_header_line, make_named_fields
from time import time
from random import random

dorky = re.compile('(\d+)([ACGT]+)>([ACGT]+)')
#########################################
def clean_cdna_change_annotation(old_annot):
    new_annot = old_annot.replace ("c.", "")
    if '>' in new_annot:
        match_return = re.match(dorky,new_annot)
        if not match_return:
            # I have no idea what this is
            new_annot = old_annot
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
def one_allele_normal(field) :
    tumor  = {field['tumor_seq_allele1'], field['tumor_seq_allele2']} # this is a set
    normal = {field['match_norm_seq_allele1'], field['match_norm_seq_allele2']}
    return len(tumor&normal)>0

#########################################
def is_useful(fields, header):
    non_info  = ['missing', '', '.', '-']
    return fields != None and fields.has_key(header) and fields[header] != None and  not fields[header].replace(" ", "") in non_info

#########################################
def new_conflict_annotation (base_annot, list_of_ids):
    annot = ""
    for db_id in list_of_ids:
        if len(annot)>0:  annot += "; "
        annot += "%s with %d" % (base_annot, db_id)
    return annot

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
def selected_info_overlap (existing_fields, new_fields, field_selection):

    # indicator (True/False) arrays
    existing_has_info = map (lambda f: existing_fields[f]is not None and len(existing_fields[f]) > 0, field_selection)
    new_has_info      = map (lambda f: new_fields[f] is not None and len(new_fields[f])>0, field_selection)
    # list of indices where both have info
    both_have_info = filter (lambda x: existing_has_info[x] and new_has_info[x], range(len(field_selection)) )
    # list of indices where both ahve info and  info is the same
    info_the_same  = filter (lambda x: existing_fields[field_selection[x]] == new_fields[field_selection[x]], both_have_info )

    # is information the same in all fields that exist in both entries?
    if len(both_have_info) != len (info_the_same):
        # nope - we have a conflict
        return "conflict"
    elif len (both_have_info) < sum (1 for x in new_has_info if x):
        # the new entry has more info
        return "new covers existing"

    return "existing covers new"

#########################################
def is_exact_duplicate (new_fields, existing_fields_by_database_id):

    exists_exact_duplicate = False
    for db_id, fields in existing_fields_by_database_id.iteritems():
        this_is_exact_duplicate = True
        for field_name, value in fields.iteritems():
            if field_name in ['id', 'conflict']: continue
            if not new_fields.has_key(field_name) or str(value) != str(new_fields[field_name]):
                this_is_exact_duplicate = False
                break
        if this_is_exact_duplicate:
            exists_exact_duplicate = True
            break
    return exists_exact_duplicate

#########################################
def diagnose_duplication_reasons (existing_fields_by_database_id, new_fields):

    diagnostics = {}
    # is this an exact duplicate by any chance?
    if is_exact_duplicate (new_fields, existing_fields_by_database_id):
        return diagnostics

    allele_fields =  ['reference_allele', 'tumor_seq_allele1', 'tumor_seq_allele2', 'match_norm_seq_allele1', 'match_norm_seq_allele2']
    interpretation_fields = ['cdna_change', 'aa_change', 'variant_classification']


    for db_id, existing_fields in existing_fields_by_database_id.iteritems():

        diagnostics[db_id] = ""
        if existing_fields['end_position'] == new_fields['end_position']:

            # giving the actual allele is the more fundamental info - go for that primarily
            allele_diagnostics = selected_info_overlap (existing_fields, new_fields, allele_fields)

            if allele_diagnostics=="conflict":
                # do we have a compound by any chance"
                if not one_allele_normal(existing_fields) and not one_allele_normal(new_fields):
                   diagnostics[db_id] = "compound heterozygous"
                else:
                   diagnostics[db_id] = "conflicting allele info"

            elif allele_diagnostics=="existing covers new":
                diagnostics[db_id] += "move on: old allele info covers"

            elif allele_diagnostics=="new covers existing":
                diagnostics[db_id] += "use new: new allele info covers"

            elif allele_diagnostics=="duplicate":
                diagnostics[db_id] += "allele info duplicate; "
                # do we have a tie breaker among the interpretation fields?
                interpretation_diagnostics = selected_info_overlap (existing_fields, new_fields, interpretation_fields)
                if interpretation_diagnostics=="conflict":
                    diagnostics[db_id] += "conflicting interpretation"
                elif interpretation_diagnostics=="existing covers new":
                    diagnostics[db_id] += "move on: old interpretation info covers"
                elif interpretation_diagnostics=="existing covers new":
                    diagnostics[db_id] += "use new: new interpretation info covers"
                elif interpretation_diagnostics=="duplicate":
                    diagnostics[db_id] += "move on: interpretation info duplicate"
                else: # we should not really be here - this is is just the default behavior so we do not crash on this
                    diagnostics[db_id] += "move on: keep old"

            else:  # we should not really be here - this is is just the default behavior so we do not crash on this
                 diagnostics[db_id] += "move on: keep old"

        else:  # the end position is not the same
            diagnostics[db_id] = "end positions different   *%s*   *%s* " % (
            existing_fields['end_position'], new_fields['end_position'])
            # how could that happen?
            # one of the possibilities (gosh, will I have to go through all of them?)
            # is that a frameshift mutation is interpreted differently
            # hard to find a robust solution
            if existing_fields['variant_classification'] == 'frame_shift_del' and one_allele_normal(
                    existing_fields) and one_allele_normal(new_fields):
                # do nothing, this is a different interpretation of the same mutation
                diagnostics[db_id] = "move on: different classification"

            elif not one_allele_normal(existing_fields) and not one_allele_normal(new_fields):
                # store, this is possibly compound  heterozygous
                diagnostics[db_id] = "compound heterozygous"
            else:
                # I don't know what this is
                diagnostics[db_id] = "conflict: different length"

    return diagnostics

 #########################################

#########################################
def store_conflicts_and_duplicates (cursor, expected_fields, new_fields):

    table = 'conflict_mutations'
    # is this an exact copy of something we have already stored as duplicate?
    sample_barcode_short = new_fields['sample_barcode_short']
    chromosome           = new_fields['chromosome']
    start_position       = new_fields['start_position']

    qry  = "select * from %s  " % table
    qry += "where sample_barcode_short = '%s' " % sample_barcode_short
    qry += "and chromosome = '%s'  " % chromosome
    qry += "and start_position = %s  " % start_position
    existing_rows = search_db(cursor, qry)

    new_entry = False
    if not existing_rows: #this is the first time we see this entry
        new_entry = True
    else:
        # take a more careful look
        existing_fields_by_database_id  = dict( zip (map (lambda x: x[0], existing_rows),  map (lambda x: make_named_fields(expected_fields, x[1:]), existing_rows) ))
        # if this is not an exact duplicate, store
        new_entry = not is_exact_duplicate(new_fields, existing_fields_by_database_id)

    if new_entry:
        insert_into_db(cursor, table, new_fields)
        return "stored in conflicts and duplicates table"
    else:
        return "found stored in conflicts and duplicates table"

#########################################
def resolve_duplicate (cursor, table, expected_fields, existing_rows, new_fields):

    existing_fields_by_database_id  = dict( zip (map (lambda x: x[0], existing_rows),  map (lambda x: make_named_fields(expected_fields, x[1:]), existing_rows) ))

    # try to diagnose how come we have multiple reports for a mutation starting at a given position
    diagnostics = diagnose_duplication_reasons (existing_fields_by_database_id, new_fields)

    # if this is an exact duplicate do nothing (this is mainly to protect us from re-runs filling the database with the same things by mistake)
    if not diagnostics: return "exact duplicate"

    if False:
        print "+"*18
        for dbid, diag in diagnostics.iteritems():
            print dbid, diag
        #   if False and new_fields['start_position'] == 35043650:
        print
        print "="*20
        for header in new_fields.keys():
            print "%30s    [ %s ] " % (header, new_fields[header]),
            for existing_fields in existing_fields_by_database_id.values():
                print " [ %s ] " % (existing_fields[header]),
            print
        print
        print

    # in the order of difficulty ...
    # 1) if it is a duplicate or subset  of one of the existing entries, we do not have to worry about it any more
    this_is_a_duplicate = len(filter(lambda x: "move on" in x, diagnostics.values() ))>0
    if this_is_a_duplicate:
        #store the new entry in the conflict_mutations table
        id_list_string = ""
        for db_id in filter(lambda x: "move on" in diagnostics[x], diagnostics.keys()):
            if len(id_list_string)>0: id_list_string+= ","
            id_list_string += str(db_id)
        descr_string = "duplicate of (" + id_list_string + ") in table %s" % table
        new_fields['conflict'] = descr_string
        return store_conflicts_and_duplicates (cursor, expected_fields, new_fields)


    there_is_no_conflict = len(filter(lambda x: "conflict" in x, diagnostics.values() ))==0

    dbids_to_be_replaced = filter(lambda dbid: "use new" in diagnostics[dbid], existing_fields_by_database_id.keys() )

    # 2) if there is no conflict, just replace the rows that the new one covers (has a superset info)
    if there_is_no_conflict:
        if len(dbids_to_be_replaced)>0:
            update_db (cursor, table, dbids_to_be_replaced[0], new_fields)
            for db_id in dbids_to_be_replaced:
                # store in the duplicates and conflicts
                store_conflicts_and_duplicates (cursor, expected_fields, existing_fields_by_database_id[db_id])
            for db_id in dbids_to_be_replaced[1:]:
                # delete from the main table
                qry = "delete from %s where id=%d" % (table, db_id)
                search_db (cursor, qry)
            return "superset"

        else: # there is one alternative possibility to covered entries: candidate compound mutations
            compound_db_ids = filter (lambda dbid: "compound" in diagnostics[dbid], existing_fields_by_database_id.keys())
            if len(compound_db_ids)>0:
                new_fields['conflict'] = "compound"
                insert_into_db (cursor, table, new_fields)
                return "compound"
            else:
                print " *** ",  filter(lambda dbid: "use new" in diagnostics[dbid], existing_fields_by_database_id.keys() )
                print "we should not have ended here ..."
                exit(1)

    # 3) if there is a conflict, but no rows to replace, add the new row and mark conflict
    #    (potentially compound get the similar treatment)
    conflicted_db_ids = filter(lambda dbid: "conflict" in diagnostics[dbid], existing_fields_by_database_id.keys() )
    if len(dbids_to_be_replaced)==0:
        new_id = insert_into_db (cursor, table, new_fields)
        for db_id in  conflicted_db_ids:
            update_conflict_field (cursor, table, db_id, existing_fields_by_database_id[db_id], new_id, "unresolved")
            update_conflict_field (cursor, table, new_id, new_fields, db_id, "unresolved")

    # 4) if there there are rows to replace and there is conflict with some other rows,
    #    replace the rows and reconsider conflicts
    else:
        # new entry
        id_to_reuse =  dbids_to_be_replaced[0]
        new_fields['conflict'] = new_conflict_annotation ("unresolved", conflicted_db_ids)
        update_db (cursor, table, id_to_reuse, new_fields)

        for db_id in dbids_to_be_replaced:
            # store in the duplicates and conflicts
            store_conflicts_and_duplicates (cursor, expected_fields, existing_fields_by_database_id[db_id])

        for db_id in dbids_to_be_replaced[1:]:
            qry = "delete from %s where id=%d" % (table, db_id)
            search_db (cursor, qry)
        # corrected old entries
        for db_id in conflicted_db_ids:
            new_annot =  new_conflict_annotation ("unresolved", map(lambda x: id_to_reuse if x==db_id else x, conflicted_db_ids))
            update_db (cursor, table, db_id, {'conflict': new_annot})

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

    #db_names  = ["BLCA"]

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
