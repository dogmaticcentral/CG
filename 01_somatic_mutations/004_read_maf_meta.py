#!/usr/bin/python -u
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
# store the meta info about the maf files: name and the reference genome,  for now

import os
from tcga_utils.mysql import *
from tcga_utils.utils import *
from random import random
import urllib2
from HTMLParser import HTMLParser
from bs4 import BeautifulSoup

##################################################################################
##################################################################################
def main():
    db = connect_to_mysql()
    cursor = db.cursor()
    db_names = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
                "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "REA",
                "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]


    for db_name in db_names:
        # check db exists
        qry = "show databases like '%s'" % db_name
        rows = search_db(cursor, qry)
        if not rows:
            print db_name, "not found"
            continue

        print " ** ", db_name
        switch_to_db(cursor, db_name)

        db_dir = "/mnt/databases/TCGA/" + db_name + "/Somatic_Mutations"
        if not os.path.isdir(db_dir):
            print "directory " + db_dir + " not found"
            exit(1)

        qry = "select * from mutations_meta"
        rows = search_db(cursor, qry)
        if not rows:
            print "no meta info found"
            continue

        maf_file = {}
        maf_diagnostics = {}
        for row in rows:
            [meta_id, file_name, quality_check, assembly, diagnostics] = row
            #if diagnostics and "tumor alleles identical" in diagnostics:
            if True:
                print "\t %4d  %50s   " % (meta_id, file_name)
                print "\t\t %6s   %6s      %s" % (quality_check, assembly, diagnostics)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
