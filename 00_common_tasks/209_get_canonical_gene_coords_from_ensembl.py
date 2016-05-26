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
# ask ensembl which genome assemblies it has available
# "hg19 is typically thought of as a subset of GRCh37", whatever that means (https://www.biostars.org/p/123343/)
# it lools like NCBI36 correponds to hg18 (https://www.biostars.org/p/17107/)
from tcga_utils.mysql import *



#########################################
def get_latest_ensembl_db_for_assembly(cursor):

    qry = "show databases like 'homo_sapiens_core_%'"
    rows = search_db(cursor, qry)
    if not rows:
        print "no response to ", qry
        exit(1)
    human_dbs = [row[0] for row in rows if not 'express' in row[0]]

    latest = {}
    for human_db in human_dbs:
        qry = " select meta_value  from %s.meta   where meta_key='assembly.default'" % human_db
        rows = search_db(cursor, qry)
        if not rows:
            print "no response to ", qry
        else:
            assembly = rows[0][0]
            latest[assembly]  = human_db

    return latest

#########################################
def main():
    db     = connect_to_mysql(conf_file="/home/ivana/.ensembl_mysql_conf")
    if not db:
        print "failed opening ensembl mysql connection"
        exit(1)
    cursor = db.cursor()
    latest = get_latest_ensembl_db_for_assembly(cursor)
    for assembly, human_db in latest.iteritems():
        print "latest db for %s is %s" % (assembly, human_db)

    cursor.close()
    db.close()

#########################################
if __name__ == '__main__':
    main()
