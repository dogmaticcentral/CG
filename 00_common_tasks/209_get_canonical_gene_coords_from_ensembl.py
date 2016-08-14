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

    db_local = connect_to_mysql()
    cursor_local = db_local.cursor()
    switch_to_db(cursor_local, 'name_resolution')

    db     = connect_to_mysql(conf_file="/home/ivana/.ensembl_mysql_conf")
    if not db:
        print "failed opening ensembl mysql connection"
        exit(1)
    cursor = db.cursor()
    latest = get_latest_ensembl_db_for_assembly(cursor)

    chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y']

    fields_to_dwld = ["hugo_name", "stable_id", "biotype", "strand", "txStart", "txEnd", "exonStarts",  "exonEnds"]

    
    for assembly, human_db in latest.iteritems():
        if assembly != 'NCBI36': continue
        print "latest db for %s is %s" % (assembly, human_db)
        target_dir = "/mnt/databases/ensembl/canonical_gene_coords/" + assembly
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)

        switch_to_db(cursor,  human_db)
        # first, lets find the chromosomes
        qry = "select * from  coord_system where name= 'chromosome' and attrib='default_version'"
        rows = search_db(cursor, qry)
        if not rows or len(rows)>1:
            print "something went wrong with ", qry
            exit(1)
              
        [coord_system_id, species_id, name, version, rank, attrib] =  rows[0]
        print "\t", coord_system_id, species_id, name, version, rank, attrib
        if species_id != 1 or version!=assembly:
            print "oink ?!"
            exit(1)

        qry = "select *  from seq_region where coord_system_id=%d" % coord_system_id
        rows = search_db(cursor, qry)
        region_id = {}
        for row in rows:
            [seq_region_id, name, coord_system_id, length] = row
            if name in chromosomes:
                region_id[name] = seq_region_id


        for chrom in chromosomes:
            
            #################################################################
            outf = open(target_dir + "/chr" + chrom+".csv", "w")
            print  >>outf,  "\t".join( fields_to_dwld )

            # gene:  gene_id, canonical_transcript_id
            # older verions keep stable id in a different table
            if 'NCBI' in assembly:
                qry  = "select distinct gene.gene_id,  gene_stable_id.stable_id, gene.canonical_transcript_id, gene.biotype "
                qry += "from gene, gene_stable_id "
                qry += "where gene.seq_region_id=%d " % region_id[chrom]
                qry += "and gene.gene_id = gene_stable_id.gene_id"
            else:
                qry = "select distinct gene_id,  stable_id, canonical_transcript_id, biotype from gene where seq_region_id=%d" % region_id[chrom]
            rows = search_db(cursor, qry)

            gene_id = {}
            canonical_transcript_id = {}
            biotype = {}

            for row in rows:
                # tanslation start and end!
                [g, stable_id, c, b] = row
                gene_id[stable_id] = g
                canonical_transcript_id[stable_id] = c
                biotype[stable_id] = b

            for stable_id in gene_id.keys():
                qry = "select symbol from hgnc where ensembl_gene_id = '%s'" % stable_id
                rows = search_db(cursor_local, qry)
                if not rows:
                    hugo_symbol = "not_found"
                else:
                    hugo_symbol = ",".join([row[0] for row in rows]) # there should not be more than one, but oh well
                this_is_protein_coding = (biotype[stable_id] == 'protein_coding')

                # first we establish in which exon does the translation start, and in which it ends
                if this_is_protein_coding:
                    # translation start and translation end
                    # seq_start = 1-based offset into the relative coordinate system of start_exon_id
                    # seq_end   = 1-based offset into the relative coordinate system of end_exon_id
                    qry  = "select seq_start, start_exon_id, seq_end, end_exon_id from translation "
                    qry += "where transcript_id = %d" %  canonical_transcript_id[stable_id]
                    rows = search_db(cursor, qry)
                    if not rows or len(rows) > 1:
                        print "something off for", stable_id, "\n"
                        print qry
                        print rows
                        exit(1)
                    [seq_start, start_exon_id, seq_end, end_exon_id] = rows[0]
                else:
                    [seq_start, start_exon_id, seq_end, end_exon_id] = [-1]*4


                # then we download all exons
                qry  = "select e.exon_id, e.seq_region_id, e.seq_region_start, e.seq_region_end, e.seq_region_strand, "
                qry += "e.phase, e.end_phase, e.is_current, e2t.rank  "
                qry += "from exon as e,  exon_transcript as e2t "
                qry += "where  e.exon_id = e2t.exon_id "
                qry += "and e2t.transcript_id = %d " %  canonical_transcript_id[stable_id]

                rows = search_db(cursor, qry)
                exon_starts = ['0'] *len(rows)
                exon_ends   = ['0'] *len(rows)
                tx_start    = -1
                tx_end      = -1

                for row in rows:
                   [exon_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, phase, end_phase, is_current, rank]  = row
                   # sanity checking or paranoia?
                   if region_id[chrom] !=  seq_region_id:
                       print "%s: seq region id mismatch (!?)   gene %d   exon %d " % (stable_id, region_id[chrom], seq_region_id)
                       exit(1)
                   exon_starts[rank-1] = str(seq_region_start)
                   exon_ends  [rank-1] = str(seq_region_end)

                   if this_is_protein_coding:
                       if exon_id == start_exon_id:
                           if seq_region_strand > 0:
                                tx_start = seq_region_start + seq_start -1
                           else:
                                tx_end   = seq_region_end  - (seq_start-1)

                       if exon_id == end_exon_id:
                           if seq_region_strand > 0:
                                tx_end = seq_region_start + seq_end-1
                           else:
                                tx_start = seq_region_end - (seq_end-1)

                if tx_end>0 and  tx_start > tx_end:
                    tx_start,tx_end = tx_end, tx_start
                # >> outf,
                print >> outf, "\t".join([hugo_symbol, stable_id, biotype[stable_id], str(seq_region_strand),
                                          str(tx_start), str(tx_end), ",".join(exon_starts), ",".join(exon_ends)])

    cursor.close()
    db.close()

    cursor_local.close()
    db_local.close()



#########################################
if __name__ == '__main__':
    main()
