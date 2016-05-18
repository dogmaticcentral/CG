#!/usr/bin/python
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
from tcga_utils.mysql import *
from tcga_utils.ucsc import *
import os

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
    for assembly in ["hg18", "hg19"]:
        switch_to_db(cursor, assembly) # mouse build name
        # no Y chromosome, we are looking at uterus tissue
        chromosomes = ["chr"+str(x) for x in range(1,23)] + ["chrX", "chrY"]
        fields_to_dwld = ["hugo_name", "ucsc_transcript_id", "ucsc_protein_id","strand", "txStart", "txEnd", "exonStarts",  "exonEnds"]
        target_dir = "/mnt/databases/ucsc/canonical_gene_coords/" + assembly

        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        for chrom in chromosomes:
            print "downloading data for", chrom
            outf = open(target_dir + "/" + chrom+".csv", "w")
            print  >>outf,  "\t".join( fields_to_dwld )
            # get canonical transcripts
            qry = "select  chromStart, chromEnd, transcript, protein from knownCanonical where chrom='%s' " % chrom
            rows = search_db(cursor, qry)
            if not rows:
                print "no return from ", qry
                exit(1)
            print "number of entries: ", len(rows)
            for row in rows:
                if len(row) < 4:
                    print row
                    exit(1)
                [tx_start, tx_end, transcript_id, protein_id] = row

                # first lets find the exons boundaries,; they are the actual definition of a splice
                qry = "select strand, txStart, txEnd, exonStarts, exonEnds from knownGene where name='%s' " % (transcript_id)
                rows2 = search_db(cursor, qry)
                if not rows2:
                    print "no return from ", qry
                    exit(1)
                if len(rows2)>1:
                    print "more than one return for the canonical transcript ?"
                    exit(1)
                [strand, tx_start, tx_end, exon_starts, exon_ends] = rows2[0]

                # now let's see if we have a name to with this splice
                qry = "select distinct name2 from refGene where exonStarts='%s' and exonEnds='%s' " % (exon_starts, exon_ends)
                qry += "and chrom = '%s' " % chrom
                rows2 = search_db(cursor, qry)
                hugo_name = "not_found"
                if rows2:
                    if len(rows2) > 1:
                        #"different HUGO name for the same gene(?)"
                        # apprently yes, esp for orfs
                        hugo_name = ",".join( [row2[0] for row2 in rows2])
                    else:
                        hugo_name = rows2[0][0]

                # lets save the transcript and protein ID to use later to get the peptide and mRNA sequence
                print >>outf, "\t".join([hugo_name, transcript_id, protein_id, strand, str(tx_start), str(tx_end), exon_starts, exon_ends])


    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


