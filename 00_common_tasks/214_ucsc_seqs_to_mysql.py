#!/usr/bin/python
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
from tcga_utils.utils import make_named_fields
from tcga_utils.ucsc import *

import os

def make_canonical_transcript_table (cursor, db_name, tbl):
    switch_to_db (cursor, db_name)

    qry = ""
    qry += "  CREATE TABLE  %s (" % tbl
    qry += "     id INT NOT NULL AUTO_INCREMENT, "
    qry += "  	 transcript_id VARCHAR (50) NOT NULL, "
    qry += "  	 protein_id VARCHAR (50), "
    qry += "  	 hugo_name VARCHAR (50), "
    qry += "	 strand VARCHAR (5) NOT NULL, "
    qry += "	 tx_start INT  NOT NULL, "
    qry += "	 tx_end INT NOT NULL, "
    qry += "	 exon_starts BLOB, "
    qry += "	 exon_ends BLOB, "
    qry += "	 mrna BLOB, "
    qry += "	 protein BLOB, "
    qry += "	 PRIMARY KEY (id) "
    qry += ") ENGINE=MyISAM"

    rows = search_db(cursor, qry)
    print qry
    print rows


#########################################
def check_dir(dir):
    if not os.path.isdir(dir):
        print dir, "not found"
        exit(1)

#########################################
def check_file(file_name):

    if not os.path.isfile(file_name):
        print file_name, "not found"
        exit(1)
    if os.stat(file_name).st_size== 0:
        print file_name, "seems to be empty"
        exit(1)

#########################################
def read_coords_file(file):

    hugo_name = {}
    trancript2protein = {}
    coords = {}

    #header_line   = file.readline()
    #header_fields = header_line.split("\t")
    #print header_fields
    for line in  file.readlines()[1:]:
        fields = line.split("\t")
        # assume the fields are
        # ['hugo_name', 'ucsc_transcript_id', 'ucsc_protein_id', 'strand', 'txStart', 'txEnd', 'exonStarts', 'exonEnds']
        ucsc_transcript_id = fields[1]
        # these should be canonical transcripts, thus one for each name
        if fields[0]== "not_found":
            hugo_name[ucsc_transcript_id] = None
        else:
            hugo_name[ucsc_transcript_id] = fields[0]
        trancript2protein[ucsc_transcript_id] = fields[2]
        coords[ucsc_transcript_id]  = fields[3:]

    return [hugo_name, trancript2protein, coords]

#########################################
def read_fasta(file):
    seq = {}
    name = ""
    for raw_line in file:
       line = raw_line.strip()
       if '>' in line:
            name = line.replace('>', '')
            seq[name] = ""
       else:
           seq[name] += line
    return seq

#########################################
def main():
    # note the skip-auto-rehash option in .ucsc_myql_conf
    # it is the equivalent to -A on the mysql command line
    # means: no autocompletion, which makes mysql get up mych faster

    db     = connect_to_mysql()
    cursor = db.cursor()

    db_name = 'ucsc'
    switch_to_db(cursor, db_name)
    chromosomes = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]
    table= {"hg18":"canonical_transcripts_hg18", "hg19":"canonical_transcripts_hg19"}
    #check that both tables exit, and if not, make them
    for tbl in  table.values():
        if ( check_table_exists (cursor, db_name, tbl)):
            print tbl, " found in ", db_name
        else:
            make_canonical_transcript_table (cursor, db_name, tbl)

    for assembly in ["hg18", "hg19"]:
        print "assembly", assembly
        coordinates_dir = "/mnt/databases/ucsc/canonical_gene_coords/" + assembly
        check_dir(coordinates_dir) # will die here if nonexistent
        target_dir = {}
        for seq_type in ['mrna', 'pep']:
            target_dir[seq_type] = "/mnt/databases/ucsc/sequences/" + assembly + "/" + seq_type
            check_dir (target_dir[seq_type])

        seqs = {}
        for chrom in chromosomes:
            print "\t chromosome", chrom
            file_name = coordinates_dir + "/" + chrom+".csv"
            check_file(file_name)
            coords_file = open(file_name, "r")

            [hugo_name, trancript2protein, coords] = read_coords_file(coords_file)
            coords_file.close()
            for seq_type in ['mrna', 'pep']:
                file_name = target_dir[seq_type] + "/" + chrom + ".fasta"
                check_file (file_name)
                seq_file  = open(file_name, "r")
                seqs[seq_type] = read_fasta(seq_file)
                seq_file.close()

            for transcript_id, protein_id in trancript2protein.iteritems():

                fixed_fields = {'transcript_id': transcript_id}

                update_fields = make_named_fields (["strand",  "tx_start", "tx_end", "exon_starts", "exon_ends"], coords[transcript_id])
                update_fields ['protein_id'] = transcript_id
                update_fields ['hugo_name']  = hugo_name[transcript_id]

                if seqs['mrna'].has_key(transcript_id):
                    update_fields ['mrna'] = seqs['mrna'][transcript_id]
                else:
                    update_fields ['mrna'] = None

                if seqs['pep'].has_key(protein_id):
                    update_fields ['protein'] = seqs['pep'][protein_id]
                else:
                    update_fields ['protein'] = None

                store_or_update (cursor, table[assembly], fixed_fields, update_fields)

    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


