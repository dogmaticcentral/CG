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

    chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y']
    
    for assembly, human_db in latest.iteritems():
        print "latest db for %s is %s" % (assembly, human_db)
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
            
            # gene:  gene_id, canonical_transcript_id
            # older verions keep stable id in a different table
            if 'NCBI' in assembly:
                qry  = "select distinct gene.gene_id,  gene_stable_id.stable_id, gene.canonical_transcript_id from gene, gene_stable_id "
                qry += "where gene.seq_region_id=%d " % region_id[chrom]
                qry += "and gene.gene_id = gene_stable_id.gene_id"
            else:
                qry = "select distinct gene_id,  stable_id, canonical_transcript_id from gene where seq_region_id=%d" % region_id[chrom]
            rows = search_db(cursor, qry)
            print "chrom ", chrom, "number of genes", len(rows)
            continue
            gene_id = {}
            canonical_transcript_id = {}
            for row in rows:
                [g, stable_id, c] = row
                gene_id[stable_id] = g
                canonical_transcript_id[stable_id] = c

            for stable_id in gene_id.keys():
                print stable_id
                qry = "select exon.*, exon_transcript.rank  from exon,  exon_transcript "
                qry += "where  exon.exon_id = exon_transcript.exon_id "
                qry += "and exon_transcript.transcript_id = %d " %  canonical_transcript_id[stable_id]
                print qry
                rows = search_db(cursor, qry)
                for row in rows:
                    print row
                exit(1)
                
        # exon_transcript: exon_id, transcript_id, rank
        # exon: exon_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, phase, end_phase, is_current
    cursor.close()
    db.close()

#########################################
if __name__ == '__main__':
    main()
