#!/usr/bin/python
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
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
from mysql import *
from tcga_utils.utils import make_named_fields
from tcga_utils.ucsc import *

import os
# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

########################################
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

    transcript2hugo   = {}
    trancript2protein = {}
    coords = {}
    orphan_hugo_symbols = []

    for line in  file.readlines()[1:]:
        fields = line.rstrip().split("\t")
        # assume the fields are
        # ['transcript2hugo', 'ucsc_transcript_id', 'ucsc_protein_id', 'strand', 'txStart', 'txEnd', 'exonStarts', 'exonEnds']
        hugo =  fields[0]
        ucsc_transcript_id = fields[1]
        if hugo =='not_found' and ucsc_transcript_id=='not_found': continue # this should not really have happened
        # these should be canonical transcripts, thus one for each name
        if ucsc_transcript_id=='not_found':
            orphan_hugo_symbols.append(hugo)
            coords[hugo] = fields[3:]
        else:
            if hugo == "not_found":
                transcript2hugo[ucsc_transcript_id] = None
            else:
                transcript2hugo[ucsc_transcript_id] = fields[0]
            trancript2protein[ucsc_transcript_id] = fields[2]
            coords[ucsc_transcript_id]  = fields[3:]

    return  [orphan_hugo_symbols, transcript2hugo, trancript2protein, coords]

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
def  handle_orphan_hugo_name(cursor, table_name, assembly, chrom,  hugo_symbol, coordinates):

    [fixed_fields, update_fields] = [ {}, {}]
    [ strand, tx_start, tx_end, exon_starts, exon_ends] = coordinates


    start = [int(x) for x in tx_start.split(";")]
    end   = [int(x) for x in tx_end.split(";")]
    e_starts = exon_starts.split(";")
    e_ends   = exon_ends.split(";")
    number_of_splices = len(start)
    if len(end) != number_of_splices or  len(e_starts) != number_of_splices  or  len(e_ends) != number_of_splices:
        print hugo_symbol, "mismatch in the number of splice coordinates"
        exit(1)

    for i in range(number_of_splices):
        fixed_fields = {'hugo_symbol': hugo_symbol}
        update_fields =  make_named_fields (["strand",  "tx_start", "tx_end", "exon_starts", "exon_ends"], \
                                            [strand, start[i], end[i], e_starts[i], e_ends[i]] )

        # get the mRNA seq directly from UCSC:
        seq = segment_from_das(assembly, chrom, start[i], end[i])
        update_fields ['mrna'] = seq
        # dash "-" is appparently a naturally occuring read-through
        if seq and not 'MIR' in hugo_symbol and not 'SNOR' in hugo_symbol and not '-' in hugo_symbol:
            # store protein sequence only if it translates cleanly
            es =  [int(x)-start[i] for x in e_starts[i].split(",")] # exons starts for this particular splice
            ee =  [int(x)-start[i]+1 for x in e_ends[i].split(",")] # upper bound on the range in python
            number_of_exons = len(es)

            if len(ee) !=  number_of_exons:
                pepseq = None  # I'm not going there
            else:
                dna = "".join ( [seq[es[n]:ee[n]] for n in range(number_of_exons)] )
                if len(dna)%3==0: # again, I am not going to deal with this here, there is just to much crap in this database
                    dnaseq = Seq (dna, generic_dna)
                    if strand=="-":
                        dnaseq = dnaseq.reverse_complement()
                    pepseq = str(dnaseq.translate())
                    if  strand=="+" :
                        print hugo_symbol, start[i], end[i], e_starts[i], e_ends[i], strand
                        print pepseq
                        print
                        #exit(1)
                else:
                    pepseq = None
            #exit(1)
        #store_or_update (cursor, table[assembly], fixed_fields, update_fields)

    return

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
            table_name = "canonical_transcripts_%s_%s" % (assembly, chrom)
            file_name  = coordinates_dir + "/" + chrom+".csv"
            check_file(file_name)
            coords_file = open(file_name, "r")

            [orphan_hugo_symbols, transcript2hugo, trancript2protein, coords] = read_coords_file(coords_file)
            coords_file.close()

            for hugo_symbol in orphan_hugo_symbols:
                handle_orphan_hugo_name(cursor, table_name, assembly, chrom, hugo_symbol, coords[hugo_symbol])

            exit(1)

            for seq_type in ['mrna', 'pep']:
                file_name = target_dir[seq_type] + "/" + chrom + ".fasta"
                check_file (file_name)
                seq_file  = open(file_name, "r")
                seqs[seq_type] = read_fasta(seq_file)
                seq_file.close()


            for transcript_id, protein_id in trancript2protein.iteritems():

                fixed_fields = {'transcript_id': transcript_id}
                update_fields =  make_named_fields (["strand",  "tx_start", "tx_end", "exon_starts", "exon_ends"], coords[transcript_id])
                update_fields ['protein_id']  = transcript_id
                update_fields ['transcript2hugo'] = transcript2hugo[transcript_id]

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


