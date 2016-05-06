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
from tcga_utils.utils  import  get_expected_fields, is_useful, make_named_fields
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

    for db_name in db_names:
        print " ** ", db_name
        switch_to_db (cursor, db_name)
        if not check_table_exists (cursor, db_name, table):
            print table, " table not found in ", db_name
            continue
        qry = "select count(1) from somatic_mutations"
        rows = search_db(cursor, qry)
        if not rows or rows[0][0] ==0 :  break
        print "somatic mutations all: ", rows[0][0],
        qry = "select count(1) from somatic_mutations where conflict is not null"
        rows = search_db(cursor, qry)
        print "conflicts: ", rows[0][0]
        print

    exit(1)

    db_names  = ["KIRC"]

    for db_name in db_names:
        # check db exists
        qry = "show databases like '%s'" % db_name
        rows = search_db(cursor, qry)
        if not rows:
            print db_name, "not found"
            exit(1) # db not found

        print " ** ", db_name
        switch_to_db (cursor, db_name)
        if not check_table_exists (cursor, db_name, table):
            print table, " table not found in ", db_name
            continue

        expected_fields = get_expected_fields(cursor, db_name, table)
        # get all conflicted groups
        qry = "select * from somatic_mutations where conflict is not null"
        rows = search_db(cursor, qry)
        if not rows:
            continue
        existing_fields_by_database_id = dict(zip(map(lambda x: int(x[0]), rows),
                                                  map(lambda x: make_named_fields(expected_fields, x[1:]), rows)))

        bags = []
        keyset = set (existing_fields_by_database_id.keys())
        for db_id, fields in existing_fields_by_database_id.iteritems():
            conflict_ids = set([int(a.split(" with ")[-1]) for a in fields['conflict'].split(";")])
            if not conflict_ids <= keyset:
                print "error - conflicting ids not reciprocally labelled or not in the database"
                exit(1)
            new_bag = set(conflict_ids) | set([db_id])  # set unioin
            new_bags = []
            for bag in bags:
                if not bag & new_bag: # intersection is empty
                    new_bags.append(bag)
                else:
                    new_bag |= bag  # add the elements of this bag to the new bag
            new_bags.append(new_bag)
            bags = new_bags

        count = {}
        count['dnp/snp']  = 0
        count['suspicious normal'] = 0
        ef = existing_fields_by_database_id # I neeed a shorthand
        for bag in bags: # bag is a collection of conflicting ids
            # now comes the random collection of reasons why this duplicate might exist:
            # 1) is this a dnp rather than snp? (I choose to believe them, otherwise I'll go crazy)
            diagnosed = False
            if len(bag)==2:
                if set(map(lambda y: ef[y]['variant_type'], [x for x in bag if ef[x].has_key('variant_type')])) == set(['snp','dnp']):
                         count['dnp/snp'] += 1
                         diagnosed = True
                elif set(is_useful(ef[db_id],'match_norm_seq_allele1')for db_id in bag) == set([True, False]):
                    count['suspicious normal'] += 1
                    diagnosed = True
                elif set(is_useful(ef[db_id],'match_norm_seq_allele2')for db_id in bag) == set([True, False]):
                    count['suspicious normal'] += 1
                    diagnosed = True

            if not diagnosed:
                print
                for field in expected_fields:
                    print field, "     ",
                    for db_id in bag:
                        print ef[db_id][field], "  ",
                    print

                print

        print " number of conflicting groups = %d" % len(bags)
        for k, v in count.iteritems():
            print k, v
        exit(1)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
