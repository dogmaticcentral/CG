#!/usr/bin/python -u
# store the meta info about the maf files: name and the reference genome,  for now
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

from old_tcga_tools.tcga_utils.utils import *
from random import random
from old_tcga_tools.tcga_utils.ucsc import segment_from_das



##################################################################################
##################################################################################
def main():
    db = connect_to_mysql()
    cursor = db.cursor()
    db_names = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
                "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "REA",
                "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    for db_name in db_names:

        switch_to_db(cursor, db_name)
        qry  = "select * from mutations_meta where diagnostics like '%stutter%'"
        rows = search_db(cursor,qry)
        if not rows: continue
        print
        print " ** ", db_name
        for row in rows:
            terms = row[-1].split(";")
            for term in terms:
                if not "stutter" in term: continue
                sample_ids = term.split("=>")[-1].replace(" ","").split(",")
                sample_ids_quoted = map[lambda x: '"' + x + '"', sample_ids].join(",")
                print sample_ids_quoted
                #qry = "delete from somatic_mutations where tumor_sample_barcode in (%s) " % sample_ids
                #search_db(cursor,qry)
                #qry = "delete from metastatic_mutations where tumor_sample_barcode in (%s) " % sample_ids
                3search_db(cursor,qry)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
