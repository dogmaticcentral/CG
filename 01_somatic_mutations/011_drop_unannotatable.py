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

# if we still have some cases that are un0annotated at the protein level - just drop them

import os.path
import re, commands
from tcga_utils.mysql import  *
from tcga_utils.utils import  get_expected_fields, is_useful, make_named_fields
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
    db_names  = ["CHOL"]
    conflicts = {}
    for db_name in db_names:
        switch_to_db (cursor, db_name)
        if not check_table_exists (cursor, db_name, table):
            print table, " table not found in ", db_name
            continue
        qry  = "select id, conflict from %s " % table
        qry += " where aa_change is null and"
        qry += " variant_classification ='missense_mutation' "
        rows = search_db(cursor, qry)
        if not rows or rows[0][0] == 0:  continue
        for row in rows:
            print row
            [id, conflict] = row
            if not conflict:
                qry = "delete from %s where id=%d" % (table, id)
                search_db(cursor, qry)
            else:
                pass
        # if conflict is null, then delete, otherwise be a bit more careful

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
