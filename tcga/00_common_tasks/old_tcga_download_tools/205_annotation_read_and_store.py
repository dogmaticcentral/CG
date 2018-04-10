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

# to speed up things make index:  (int 01_maf_somatic_mutations_table.sql)
# create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position)
# note: I am assuming then that no tumor will have the exact same mutation in both alleles, for the simple reason that I do nto see
# how would the two entries in the database then be distinguished from one another
# (rather if tumor_seq_allele1 == tumor_seq_allele2 != match_norm_seq_allele1  and tumor_seq_allele2 != match_norm_seq_allele2
# then I have somehting like that)


import MySQLdb
from   tcga_utils.mysql import  *
import commands


#########################################
def make_annotation_table(cursor, fields):

    qry = "";
    qry += "  CREATE TABLE annotation ("

    for field in fields:
        if field == 'id':
            qry += "  	 %s INT NOT NULL PRIMARY KEY" % field
        elif field in ['disease',  'item_barcode',]:
            qry += "  	 %s VARCHAR (50) NOT NULL" % field
        else:
            qry += "  	 %s BLOB " % field
        if not field == fields[-1]:
            qry += ", "
    qry += ") ENGINE=MyISAM"
    rows = search_db(cursor, qry)
    print qry
    print rows

    qry = "";
    qry += "create index barcode_idx on annotation (item_barcode)"
    rows = search_db(cursor, qry)
    print qry
    print rows


#########################################
def store (cursor, header_fields, fields):

    fixed_fields  = {}
    update_fields = {}

    if (len(fields) != len(header_fields)): return
    for i in range( len(header_fields) ):
        field = fields[i]
        header = header_fields[i]
        if (header in ['id'] ):
            fixed_fields[header] = int(field)
        else:
            update_fields[header] = field

    ok = store_or_update (cursor, 'annotation', fixed_fields, update_fields)
    if not ok:
        print 'store failure:'
        print fixed_fields
        print update_fields
        store_or_update (cursor, 'annotation', fixed_fields, update_fields, verbose=True)
        exit(1)



#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_name = 'tcga_meta'
    table = 'annotation'

    switch_to_db (cursor, db_name)

    table_exists = check_table_exists (cursor, db_name,  table)
    if table_exists:
        print table, " found in ", db_name
    else:
        print table, " not found in ", db_name

    annotation_file = '/mnt/databases/TCGA/tcga_annotations.txt'
    inf = open(annotation_file)
    for line in inf:
        if line[0] == '%': continue
        line = line.rstrip().replace("'", "");

        if line[:2] == 'ID':
            header_fields = [ '_'.join(x.lower().split()) for x in line.split('\t')]
            print header_fields
            if not table_exists: make_annotation_table(cursor, header_fields)
        else:
            fields = line.split('\t')
            store (cursor, header_fields, fields)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
