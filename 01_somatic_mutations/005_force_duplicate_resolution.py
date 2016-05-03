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

    db_names  = ["ACC"]

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

        expected_fields = get_expected_fields(cursor, db_name, table)
        # get all conflicted groups
        qry = "select * from somatic_mutations where conflict is not null"
        rows = search_db(cursor, qry)
        if not rows:
            continue
        existing_fields_by_database_id = dict(zip(map(lambda x: int(x[0]), rows),
                                                  map(lambda x: make_named_fields(expected_fields, x[1:]), rows)))

        bags = []
        for db_id, fields in existing_fields_by_database_id.iteritems():
            conflict_ids = [int(a.split(" with ")[-1]) for a in fields['conflict'].split(";")]
            new_bag = set(conflict_ids) | set([db_id]) # set unioin
            new_bags = []
            for bag in bags:
                if not bag & new_bag: # intersection is empty
                    new_bags.append(bag)
                else:
                    new_bag |= bag  # add the elements of this bag to the new bag
            new_bags.append(new_bag)
            bags = new_bags


        for bag in bags: # bag is a collection of conflicting ids
            print
            for field in expected_fields:
                print field, "     ",
                for db_id in bag:
                    print existing_fields_by_database_id[db_id][field], "  ",
                print
            print
        print

        print " number of conflicting groups = %d" % len(bags)
        exit(1)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
