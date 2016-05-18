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
from tcga_utils.utils  import  make_named_fields, is_informative
from tcga_utils.ucsc  import  *
from time import time

verbose = True
#########################################
def store_fixed_row (cursor, fixed_row):
    return


mutation_annot_pattern = re.compile('(\D+)(\-*\d+)(\D+)')
#########################################
def parse_mutation (mutation):

    match_return = re.match(mutation_annot_pattern, mutation)
    mut_from = match_return.group(1)
    mut_to   = match_return.group(3)
    mut_position = int (match_return.group(2))
    return [mut_position, mut_from, mut_to]

#########################################
def check_aa_type (ucsc_cursor, assembly_dict, fields):

    checks = True
    fixed_row = {}
    conflict  = fields['conflict']
    aa_change = fields['aa_change']
    variant_classification =  fields['variant_classification']
    # I'll fix the absolute minimum that I can scrape by with
    if not conflict and (variant_classification!="missense_mutation" or is_informative(aa_change)):
        return [checks, fixed_row]
    id = fields['id']
    start_position = fields['start_position']
    end_position   = fields['end_position']
    tumor1      = fields['tumor_seq_allele1']
    tumor2      = fields['tumor_seq_allele2']
    norm1       = fields['match_norm_seq_allele1']
    norm2       = fields['match_norm_seq_allele2']
    reference   = fields['reference_allele']
    aa_change   = fields['aa_change']
    cdna_change = fields['cdna_change']
    meta_info_index = fields['meta_info_index']
    assembly   = assembly_dict[meta_info_index]
    chromosome = fields['chromosome']

    ucsd_segment = segment_from_das(assembly, chromosome, start_position, end_position)
    print id, assembly, chromosome, start_position, end_position, ucsd_segment
    print reference, norm1, norm2, tumor1, tumor2, cdna_change, aa_change
    print parse_mutation (cdna_change)
    print parse_mutation (aa_change)
    print "conflict: ", conflict
    exit(1)
    return [checks, fixed_row]

#########################################
def get_assemblies (cursor):
    assembly = {}
    qry = "select id, assembly from mutations_meta"
    rows = search_db(cursor, qry)
    if not rows:
        print "assembly not found"
        exit(1) # db not found
    for row in rows:
        assembly[row[0]] = row[1]
    return assembly

#########################################
def main():

    db     = connect_to_mysql()
    cursor = db.cursor()

    ucsc = connect_to_mysql(conf_file="/home/ivana/.ucsc_mysql_conf")
    if not ucsc:
        print "failed opening ucsc mysql connection"
        exit(1)
    ucsc_cursor = db.cursor()

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

    db_names = ["LUAD"]

    chunk = 10 # we process rows 10 by 10+
    offset = -chunk
    for db_name in db_names:
        qry = "show databases like '%s'" % db_name
        rows = search_db(cursor, qry)
        if not rows:
            print db_name, "not found"
            exit(1) # db not found

        print " ** ", db_name
        switch_to_db (cursor, db_name)
        if ( check_table_exists (cursor, db_name, table)):
            print table, "table found in ", db_name
        else:
            print table, "table not found in ", db_name
        header_fields = get_column_names (cursor, db_name, table)
        if not header_fields:
            print "\t no columnn names (?)"
            continue

        assembly = get_assemblies (cursor)

        done = False
        while not done:
            offset += chunk
            if offset and not offset%1000: print "offset: ", offset
            qry = "select * from %s limit %d, %d" % (table, offset, chunk)
            rows = search_db(cursor, qry)
            if not rows:
                done = True
                continue
            for row in rows:
                [checks, fixed_row] = check_aa_type (ucsc_cursor, assembly, make_named_fields (header_fields, row) )
                if checks: continue
                store_fixed_row (cursor, fixed_row)
        break
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
