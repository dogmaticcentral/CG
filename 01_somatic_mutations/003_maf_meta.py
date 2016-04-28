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

import os
from tcga_utils.mysql import *
from tcga_utils.utils import *
from random import random
import urllib2
from HTMLParser import HTMLParser
from bs4 import BeautifulSoup


#########################################
def find_reference_genome(maffile):
    ref_gen = ""
    header_fields = process_header_line(maffile)

    ch_idx     = header_fields.index('chromosome')
    start_idx  = header_fields.index('start_position')
    end_idx    = header_fields.index('end_position')
    ref_al_idx = header_fields.index('reference_allele')

    assemblies = ["hg38", "hg19", "hg18", "hg17"]
    number_of_correct_matches = {}
    for assembly in assemblies:
        number_of_correct_matches[assembly] = 0
    # find random examples from each chromosome
    inff = open(maffile, "r")
    sample_size = 0
    line_ct = 0
    for line in inff:
        if line.isspace(): continue
        if line[0] == '#': continue
        line_ct += 1
        if line_ct <= 1: continue
        if random() > 0.01: continue
        field = line.split("\t")
        chrom = field[ch_idx]
        start = field[start_idx]
        end = field[end_idx]
        ref = field[ref_al_idx]
        # print chrom, start, end, al1, al2, norm1, norm2
        sample_size += 1
        for assembly in assemblies:
            das_request = "http://genome.ucsc.edu/cgi-bin/das/%s/" % assembly
            das_request += "dna?segment=chr%s:%s,%s" % (chrom, start, end)
            response = urllib2.urlopen(das_request)
            html = response.read()
            soup = BeautifulSoup(html, 'html.parser')
            if not soup or not soup.dna or not soup.dna.string: continue
            # print " ** ", ref.upper(), soup.dna.string.strip().upper()
            if ref.upper() == soup.dna.string.strip().upper():
                number_of_correct_matches[assembly] += 1
        if sample_size >= 100: break

    inff.close()

    max_matches = max(number_of_correct_matches.values())
    max_assemblies = filter(lambda x: number_of_correct_matches[x] == max_matches, assemblies)
    if len(max_assemblies) > 1:
        return ["fail", "multiple genome  matches:" + " ".join(max_assemblies)]

    ref_gen = max_assemblies[0]
    # try matching using ucsf server
    # use the assembly with the smallest amount of failures
    # if both assemblies fail in more htan 10% of cases - abort
    return ["pass", ref_gen]


#########################################
def check_headers(maffile, required_fields, expected_fields):
    # required_fields are the absolute minimum we need to
    # reconstruct the mutation - that they are missing  should not happen at all       
    header_fields = process_header_line(maffile)
    if len(header_fields) == 0:
        return ["fail", "no header found"]
    missing_fields = filter(lambda x: not x in header_fields, required_fields)
    if len(missing_fields) > 0:
        return ["fail", "missing required fields:  " + " ".join(missing_fields)]

    # we know we eill add tthese two filds to the original maf: 'sample_barcode_short', 'meta_info_index'
    missing_fields = filter(lambda x: not x in header_fields+['sample_barcode_short', 'meta_info_index'], expected_fields)
    if len(missing_fields) > 0:
        return ["warn", "missing expected fields:  " + " ".join(missing_fields)]

    return ["pass", ""]


##################################################################################
# checking for the following, as seen in
# broad.mit.edu_LIHC.IlluminaGA_DNASeq_automated.Level_2.1.0.0/
# An_TCGA_LIHC_External_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.maf 
#           273933        RPL5       Frame_Shift_Del         p.K270fs 
#           273933        RPL5       Frame_Shift_Del         p.K277fs 
#           273933        RPL5       Frame_Shift_Del         p.R279fs 
#           273933        RPL5       Frame_Shift_Del         p.Q282fs 
##################################################################################
def check_health(maffile):
    inff = open(maffile, "r")
    missing_fields = []
    header_fields = process_header_line(maffile)
    if len(header_fields) == 0:
        # though we shold have discovered this previously ....
        return ["fail", "no header found"]

    variantclass_index = header_fields.index('variant_classification')
    startpos_index     = header_fields.index('start_position')
    sample_barcode_idx = header_fields.index('tumor_sample_barcode')
    start_posns = {}
    for line in inff:
        if not 'Frame_Shift_Del' in line: continue
        field = line.split("\t");
        if not field[variantclass_index] == 'Frame_Shift_Del': continue
        sample_barcode = field[sample_barcode_idx]
        if not start_posns.has_key(sample_barcode): start_posns[sample_barcode] = []
        start_posns[sample_barcode].append(int(field[startpos_index]))
    inff.close()

    number_of_samples_with_stutter_count = {}
    for sample_barcode in start_posns.keys():
        prev = None
        count = 0
        for sp in start_posns[sample_barcode]:
            if prev and sp - prev < 5:
                count += 1
            prev = sp
        if not count in number_of_samples_with_stutter_count.keys(): number_of_samples_with_stutter_count[count] = 0
        number_of_samples_with_stutter_count[count] += 1

    # for stutter_count in sorted( number_of_samples_with_stutter_count.keys()):
    #    print "\t", stutter_count, number_of_samples_with_stutter_count [stutter_count]

    # I am not really sure what to do with this, so for now I will
    # only store it with a warning
    bad_counts = filter(lambda x: x > 100, number_of_samples_with_stutter_count.keys())
    if len(bad_counts) > 0:
        return ["warn", "%d samples have more that 100 stutter frameshift points" % len(bad_counts)]

    return ["pass", ""]


#########################################
def store_meta_info(cursor, bare_filename, overall_diagnostics):
    print
    print "storing meta info for ", bare_filename
    for diag in overall_diagnostics:
        print diag

    fixed_fields = {}
    update_fields = {}

    fixed_fields['file_name'] = bare_filename;
    if overall_diagnostics[-1][0] == "pass":
        update_fields['quality_check'] = "pass"
        update_fields['assembly'] = overall_diagnostics[-1][1]
        overall_diagnostics.pop()
    else:
        update_fields['quality_check'] = "fail"

    if len(overall_diagnostics) > 0:
        update_fields['diagnostics'] = ";".join(map(lambda x: ":".join(x), overall_diagnostics))

    store_or_update(cursor, "mutations_meta", fixed_fields, update_fields)
    return


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
        expected_fields = get_expected_fields(cursor, db_name, "somatic_mutations")

        for maffile in maf_files:
            bare_filename = "/".join ( maffile.split('/')[-2:])
            overall_diagnostics = []

            # first make sure that the file is not empty - in which case
            # we might have to go and check what's wiht the download
            if os.path.getsize(maffile) == 0:
                overall_diagnostics.append(["fail", "file empty"])
                ref_genome = ""
                store_meta_info(cursor, bare_filename, overall_diagnostics)
                continue

            # might come handy: is this data curated or not?
            if "automated" in bare_filename.lower():
                overall_diagnostics.append(["note", "automated"])
            elif "curated" in bare_filename.lower():
                overall_diagnostics.append(["note", "curated"])

            # check if the file contains  all the info we need and hope to have
            diagnostics = check_headers(maffile, required_fields, expected_fields)
            if diagnostics[0] != "pass":
                overall_diagnostics.append(diagnostics)
            if diagnostics[0] == "fail":
                ref_genome = ""
                store_meta_info(cursor, bare_filename, overall_diagnostics)
                continue

            # I am aware of one way in which the file can be corrupt, so I am checking for it
            diagnostics = check_health(maffile)
            if diagnostics[0] != "pass":
                overall_diagnostics.append(diagnostics)
            if diagnostics[0] == "fail":
                ref_genome = ""
                store(cursor, bare_filename, overall_diagnostics)
                continue

            diagnostics = find_reference_genome(maffile)

            overall_diagnostics.append(diagnostics)

            store_meta_info(cursor, bare_filename, overall_diagnostics)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
