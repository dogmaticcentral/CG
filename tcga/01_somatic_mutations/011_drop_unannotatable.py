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

    qry  = "select table_name from information_schema.tables "
    qry += "where table_schema='tcga' and table_name like '%_somatic_mutations'"
    tables = [field[0] for field in search_db(cursor,qry)]

    switch_to_db (cursor, "tcga")

    for table in tables:
        qry  = "select id, conflict from %s " % table
        qry += " where aa_change is null and"
        qry += " variant_classification ='missense_mutation' "
        rows = search_db(cursor, qry)
        if not rows or rows[0][0] == 0:  continue

        for row in rows:
            print row
            [id, conflict] = row
            # delete in any case
            qry = "delete from %s where id=%d" % (table, id)
            search_db(cursor, qry)
            if not conflict: continue

            # the conflicting mutation is not so conflicting any more
            conflicting_id = int(conflict.split(' ')[-1])
            qry = "select conflict from %s where id=%d" % (table, conflicting_id)
            rows2 = search_db(cursor, qry)
            if not rows2:
                print "the row with id %d does not exist" % conflicting_id
                qry = "delete from %s where id=%d" % (table, id)
                search_db(cursor, qry)
            else:
                for row2 in rows2:
                    if not row2[0]: continue
                    conflicts = row2[0].split(';')
                    new_conflicts = []
                    for confl in conflicts:
                        id2 = int(confl.split(' ')[-1])
                        if id2==id: continue
                        new_conflicts.append(confl)
                    if len(new_conflicts)==0:
                        qry = "update %s set conflict=NULL where id=%d" % (table, conflicting_id)
                    else:
                        newconfl = "; ".join(new_conflicts)
                        qry = "update %s set conflict=%s where id=%d" % (table, newconfl, conflicting_id)
                    search_db(cursor, qry)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
