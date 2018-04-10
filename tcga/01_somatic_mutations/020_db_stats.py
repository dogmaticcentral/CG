#!/usr/bin/python -u
# the fields in the tumor sample barcode are
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
# the fields in the tumor sample barcode are
# project - tissue source site (TSS)  - participant -
# source.vial - portion.analyte  - plate - (sequencing or characterization center)
#
#
# sample codes:
# 01	Primary solid Tumor	TP
# 02	Recurrent Solid Tumor	TR
# 03	Primary Blood Derived Cancer - Peripheral Blood	TB
# 04	Recurrent Blood Derived Cancer - Bone Marrow	TRBM
# 05	Additional - New Primary	TAP
# 06	Metastatic	TM
# 07	Additional Metastatic	TAM
# 08	Human Tumor Original Cells	THOC
# 09	Primary Blood Derived Cancer - Bone Marrow	TBM
# 10	Blood Derived Normal	NB
# 11	Solid Tissue Normal	NT
# 12	Buccal Cell Normal	NBC
# 13	EBV Immortalized Normal	NEBV
# 14	Bone Marrow Normal	NBM
# 20	Control Analyte	CELLC
# 40	Recurrent Blood Derived Cancer - Peripheral Blood	TRB
# 50	Cell Lines	CELL
# 60	Primary Xenograft Tissue	XP
# 61	Cell Line Derived Xenograft Tissue	XCL


import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *

verbose = False

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    full_name = read_cancer_names ()
    db_names = ["ACC", "BLCA", "BRCA", "CESC", "CHOL",  "COAD", "DLBC", "ESCA",
                "GBM", "HNSC", "KICH" ,"KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD",
                "LUSC",  "MESO", "OV",   "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    total_primary = 0
    total_metastatic = 0
    for db_name in db_names:

        switch_to_db (cursor, db_name)
        ############################
        qry = "select count(distinct tumor_sample_barcode) from somatic_mutations "
        rows = search_db(cursor, qry)
        primary_samples = int(rows[0][0])
        total_primary += primary_samples

        qry = "select count(distinct tumor_sample_barcode) from metastatic_mutations "
        rows = search_db(cursor, qry)
        metastatic_samples = int(rows[0][0])
        total_metastatic += metastatic_samples

        print db_name, full_name[db_name]
        print "primary tumor samples: %d   metastatic: %d     total: %d" % \
                  (primary_samples, metastatic_samples, primary_samples+metastatic_samples)
        print

    print "total tumor types:         ", len(db_names)
    print "      primary samples:     ", total_primary
    print "      metastatic samples:  ", total_metastatic
    print "      all samples:         ", total_metastatic+total_primary

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

