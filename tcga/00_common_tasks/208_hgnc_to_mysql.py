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

from tcga_utils.utils import *

#########################################
def make_hgnc_table (cursor, header_fields):
    # construct the mysql command to create the appropriate table
    # if id, varchar 30, otherwise blob - I don'r want to waste to much time on this
    qry = ""
    qry += "  CREATE TABLE hgnc ("
    qry += "     id INT NOT NULL AUTO_INCREMENT, "
    for field in header_fields:
        if "date_" in field: # these belong to the original dump and are not really informative
            pass
        elif field in ["ccds_id","lsdb","pubmed_id","mgd_id"]:
            qry += " %s BLOB, " %field
        elif "_id" == field[-3:]: # careful here bcs we may have _ids, which is better served as blob
            qry += " %s VARCHAR (30), " % field
        elif "db" == field[-2:]:
            qry += " %s VARCHAR (30), " % field
        else:
            qry += " %s BLOB, " %field
    qry += " PRIMARY KEY (id) "
    qry += ") ENGINE=MyISAM"
    rows = search_db(cursor, qry)
    print qry
    print rows


#########################################
def main():

    hgnc_raw_file = "/mnt/databases/hgnc_complete_set.txt"
    if not os.path.isfile(hgnc_raw_file):
        print hgnc_raw_file, " not found, or not a readable file"
        exit(1)
    if os.stat(hgnc_raw_file).st_size== 0:
        print hgnc_raw_file, " seems to be empty"
        exit(1)

    db     = connect_to_mysql()
    cursor = db.cursor()
    db_name = 'name_resolution'
    switch_to_db(cursor, db_name)


    hgnc_fh = open(hgnc_raw_file, "r")
    # we expect the first line to be a header line
    header_fields = [x.replace(".", "_").replace("-", "_") for x in hgnc_fh.readline().rstrip().split("\t")]

    if ( check_table_exists (cursor, db_name, 'hgnc')):
        print " hgnc found in ", db_name
        expected_fields = get_expected_fields(cursor, db_name, "hgnc")
        #qry = "drop table hgnc"
        #search_db(cursor, qry)
    else:
        make_hgnc_table (cursor, header_fields)
    #exit(1)
    for line in hgnc_fh.readlines()[1:]:
        fields = line.rstrip().split("\t")
        # the right hand fields might be non-exsistent
        for i in range(len(fields), len(header_fields)): fields.append("")
        named_fields = make_named_fields (header_fields, fields, expected_fields)
        fixed_fields = {"hgnc_id": named_fields["hgnc_id"]}
        del named_fields["hgnc_id"]
        update_fields = named_fields
        if not store_or_update(cursor, "hgnc", fixed_fields, update_fields):
            exit(1)
    hgnc_fh.close()

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()


