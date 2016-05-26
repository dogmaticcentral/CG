#!/usr/bin/python
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
from tcga_utils.mysql import *
from tcga_utils.ucsc import *
import os

#########################################
def find_overlapping_intervals (coord_set, hugo_coordinates, ucscs_strand, hugo_strand):
    intersecting_intervals = {}
    for hugo_name, hugo_coords in hugo_coordinates.iteritems():
        if ucscs_strand != hugo_strand[hugo_name]: continue
        for coordinates in hugo_coords:
            hg  = list(coordinates)
            for ucscs_coords in coord_set:
                cs = list(ucscs_coords)
                intersection = not ((cs[1] < hg[0]) or (hg[1] < cs[0]) )
                if intersection: intersecting_intervals[hugo_name] = coordinates
    return intersecting_intervals

#########################################
def check_rows (rows, qry):
    if not rows:
        print "no return from ", qry
        exit(1)
def exon_list_cleanup(exon_list):
    el = exon_list.strip()
    if el[-1]==",": el = el[:-1]
    return el

#########################################
def main():
    # note the skip-auto-rehash option in .ucsc_myql_conf
    # it is the equivalent to -A on the mysql command line
    # means: no autocompletion, which makes mysql get up mych faster

    db     = connect_to_mysql(conf_file="/home/ivana/.ucsc_mysql_conf")
    if not db:
        print "failed opening ucsc mysql connection"
        exit(1)
    cursor = db.cursor()

    chromosomes = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]
    for assembly in ["hg18", "hg19", "hg38"]:
        switch_to_db(cursor, assembly) # mouse build name

        fields_to_dwld = ["hugo_name", "ucsc_transcript_id", "ucsc_protein_id","strand", "txStart", "txEnd", "exonStarts",  "exonEnds"]

        target_dir = "/mnt/databases/ucsc/canonical_gene_coords/" + assembly
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)

        for chrom in chromosomes:
            print "downloading data for", chrom

            # get all refSeq entries for the crhomosome - they conatina hugo names
            qry = "select name2, strand, txStart, txEnd from refGene where chrom = '%s' " % chrom
            rows = search_db(cursor, qry)
            check_rows (rows, qry)

            hugo_coordinates = {}
            hugo_strand = {}
            for row in rows:
                [name2, strand, txStart, txEnd] = row
                name2 = name2.strip()
                hugo_strand[name2] = strand
                if not hugo_coordinates.has_key(name2): hugo_coordinates[name2] = set([])
                if (txStart, txEnd) not in hugo_coordinates[name2]:
                    hugo_coordinates[name2].add((txStart, txEnd))

            # all canonical transcripts
            qry  = "select transcript, protein, chromStart, chromEnd from knownCanonical where chrom='%s' " % chrom
            rows = search_db (cursor, qry)
            check_rows(rows, qry)
            canonical_coordinates = {}
            transcript2protein = {}
            for row in rows:
                [transcript_id, protein_id, txStart, txEnd] = row
                if not canonical_coordinates.has_key(transcript_id): canonical_coordinates[transcript_id] = set([])
                if (txStart, txEnd) not in canonical_coordinates[transcript_id]:
                    canonical_coordinates[transcript_id].add((txStart, txEnd))
                    transcript2protein[transcript_id] = protein_id

            # strands for each transcript
            qry  = "select strand, name from knownGene"
            rows = search_db(cursor, qry)
            check_rows (rows, qry)
            ucsc_strand = {}
            for row in rows:
                ucsc_strand[row[1]] = row[0]

            # find exact match or common_intervals
            hugo2canonical = {}
            canonical2hugo = {}
            for transcript_id  in canonical_coordinates.keys():
                canonical2hugo[transcript_id] = []
            for hugo_name  in  hugo_coordinates.keys():
                hugo2canonical[hugo_name] = []

            for transcript_id, coord_set in canonical_coordinates.iteritems():
                # try the exact match first
                common_intervals = None
                for hugo_name, hugo_coords in hugo_coordinates.iteritems():
                    if hugo_strand[hugo_name] != ucsc_strand[transcript_id]: continue
                    if hugo_coords & coord_set:
                        common_intervals = hugo_coords & coord_set
                        canonical2hugo [transcript_id].append(hugo_name)
                        hugo2canonical [hugo_name].append(transcript_id)
                if common_intervals: continue
                # if that didn't work, try the overlap
                intersecting_intervals = find_overlapping_intervals (coord_set, hugo_coordinates, ucsc_strand[transcript_id], hugo_strand)
                for hugo_name in intersecting_intervals.keys():
                    canonical2hugo [transcript_id].append(hugo_name)
                    hugo2canonical [hugo_name].append(transcript_id)

            #################################################################
            outf = open(target_dir + "/" + chrom+".csv", "w")
            print  >>outf,  "\t".join( fields_to_dwld )

            # for hugo symbols that do not have the canonical description, use the exons from the refGene table
            orphan_hugos = filter (lambda x: hugo2canonical[x]==[],   hugo2canonical.keys() )

            qry = "select name2, txStart, txEnd, exonStarts, exonEnds from refGene where chrom = '%s' " % chrom
            rows2 = search_db(cursor, qry)
            check_rows (rows2, qry)

            orphan_coordinates = {}
            for row in rows2:
                [hugo_name, tx_st, tx_end, ex_st, ex_end] = row
                if not hugo_name in orphan_hugos: continue
                if not orphan_coordinates.has_key(hugo_name):  orphan_coordinates[hugo_name] = []
                orphan_coordinates[hugo_name].append ([tx_st, tx_end, exon_list_cleanup(ex_st), exon_list_cleanup(ex_end)])

            transcript_id = "not_found"
            protein_id = "not_found"
            for hugo_name in orphan_hugos:
                tx_starts   = []
                tx_ends     = []
                exon_starts = []
                exon_ends   = []
                for [tx_start, tx_end, ex_st, ex_end] in orphan_coordinates[hugo_name]:
                    tx_starts.append(str(tx_start))
                    tx_ends.append(str(tx_end))
                    exon_starts.append (ex_st)
                    exon_ends.append (ex_end)
                ts = ";".join(tx_starts)
                te = ";".join(tx_ends)
                es = ";".join(exon_starts)
                ee = ";".join(exon_ends)
                print >>outf, "\t".join([hugo_name, transcript_id, protein_id, hugo_strand[hugo_name], ts, te, es, ee])
            # otherwise use exons from the canonical table
            qry = "select name, exonStarts, exonEnds from knownGene where chrom='%s' " % (chrom)
            rows2 = search_db(cursor, qry)
            check_rows (rows2, qry)
            exon_starts = {}
            exon_ends = {}
            for row in rows2:
                [transcript_id, ex_st, ex_end] = row
                exon_starts[transcript_id] = exon_list_cleanup(ex_st)
                exon_ends[transcript_id]   = exon_list_cleanup(ex_end)

            for transcript_id, coord_set in canonical_coordinates.iteritems():
                if len(coord_set) != 1:
                    print transcript2protein, len(coord_set)
                    continue
                coords = list(coord_set)
                tx_start = coords[0][0]
                tx_end   = coords[0][1]
                if canonical2hugo.has_key(transcript_id) and canonical2hugo[transcript_id]:
                    hugo_names = ",".join(canonical2hugo[transcript_id])
                else:
                    hugo_names = "not_found"
                protein_id = transcript2protein[transcript_id]
                # lets save the transcript and protein ID to use later to get the peptide and mRNA sequence
                print >>outf, "\t".join([hugo_names, transcript_id, protein_id, strand, \
                                         str(tx_start), str(tx_end), exon_starts[transcript_id], exon_ends[transcript_id]])


    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


