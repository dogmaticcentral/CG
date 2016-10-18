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

        maf_files = []
        for subdir in os.listdir(db_dir):
            subfull = "/".join([db_dir, subdir])
            # we'll not deal with mitochondiral mutations for now
            mafs = filter(lambda x: x.endswith(".maf") and not 'mitochondria' in x, os.listdir(subfull))
            if len(mafs) > 1:
                print "unexpected dir structure:"
                print subfull
                print mafs
                exit(1)
            if len(mafs) == 0: continue
            maf_files.append ( "/".join([db_dir, subdir,mafs[0]]) )

        required_fields = get_required_fields()
        expected_fields = filter (lambda x: x not in ['sample_barcode_short', 'meta_info_index', 'conflict'], get_expected_fields(cursor, db_name, "somatic_mutations"))

        for maffile in maf_files:

            bare_filename = "/".join ( maffile.split('/')[-2:])
            overall_diagnostics = []

            # first make sure that the file is not empty - in which case
            # we might have to go and check what's wiht the download
            if os.path.getsize(maffile) == 0:
                overall_diagnostics.append(["fail", "file empty"])
                store_meta_info(cursor, bare_filename, overall_diagnostics)
                continue

            # check if the file contains  all the info we need and hope to have
            diagnostics = check_headers(maffile, required_fields, expected_fields)
            if diagnostics[0] != "pass":
                overall_diagnostics.append(diagnostics)
            if diagnostics[0] == "fail":
                store_meta_info(cursor, bare_filename, overall_diagnostics)
                continue

            # this one should be pass or fail only
            assembly_diag = find_reference_genome(cursor, maffile, bare_filename)
            # if we fail here,  it means it is not clear what assembly was used
            if assembly_diag[0] == "fail":
                overall_diagnostics.append(assembly_diag)
                store_meta_info(cursor, bare_filename, overall_diagnostics)
                continue

            # might come handy: is this data curated or not?
            if "automated" in bare_filename.lower():
                overall_diagnostics.append(["note", "automated"])
            elif "curated" in bare_filename.lower():
                overall_diagnostics.append(["note", "curated"])


            # I am aware of one way in which the file can be corrupt (stutter), so I am checking for it
            diagnostics = check_health(maffile)
            if diagnostics[0] != "pass":
                overall_diagnostics.append(diagnostics)

            # some seq centers fill in one of the normal alleles as '-' --> actually, it looks like it's Baylor
            # (all seq centers seem not to bother checking if the other allele is different)
            diagnostics = check_norm_allele(maffile, 1)
            if diagnostics[0] != "pass":
                 overall_diagnostics.append(diagnostics)

            diagnostics = check_norm_allele(maffile, 2)
            if diagnostics[0] != "pass":
                 overall_diagnostics.append(diagnostics)

            diagnostics = check_norm_allele(maffile, 0)
            if diagnostics[0] != "pass":
                 overall_diagnostics.append(diagnostics)

            diagnostics = check_tumor_alleles(maffile)
            if diagnostics[0] != "pass":
                 overall_diagnostics.append(diagnostics)

            overall_diagnostics.append(diagnostics)

            # as a convention, (or convenience) assembly info last
            overall_diagnostics.append(assembly_diag)

            store_meta_info(cursor, bare_filename, overall_diagnostics)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
