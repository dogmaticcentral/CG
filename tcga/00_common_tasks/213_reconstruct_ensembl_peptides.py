#!/usr/bin/python

from tcga_utils.ucsc import *

import os
# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

verbose = False

#########################################
def main():
    grc = {'hg38':'GRCh38', 'hg19':'GRCh37', 'hg18':'NCBI36'}
    chromosomes = ["chr"+str(x) for x in range(1,23)] + ["chrX", "chrY"]
    for assembly in ["hg18","hg19"]:
        source_dir = "/mnt/databases/ensembl/canonical_gene_coords/" + grc[assembly]
        target_dir = {}
        for seq_type in ['mrna', 'pep']:
            target_dir[seq_type] = "/mnt/databases/ensembl/sequences/" + assembly + "/" + seq_type
            if not os.path.isdir(target_dir[seq_type]):
                os.makedirs(target_dir[seq_type])
        for chrom in chromosomes:
            print assembly, chrom
            inf =  open(source_dir+ "/" + chrom+".csv", "r")
            outf_mrna = open(target_dir['mrna'] + "/" + chrom + ".fasta", "w")
            outf_pep  = open(target_dir['pep'] + "/" + chrom + ".fasta", "w")
            for line in inf.readlines()[1:]:
                fields = line.rstrip().split( "\t")

                [hugo, stable_id, biotype, strand, tx_start, tx_end, exon_starts, exon_ends] = fields
                strand = int(strand)
                tx_start = int(tx_start)
                tx_end   = int (tx_end)

                es =  sorted([int(x) for x in exon_starts.split(",")] )
                ee =  sorted([int(x) for x in exon_ends.split(",")] )
                min_pos = es[0]
                max_pos = ee[-1]
                seq = segment_from_das(assembly, chrom, min_pos, max_pos, verbose=verbose)

                #now, let's move the coord_system
                es = [x-min_pos for x in es]
                ee = [x-min_pos for x in ee]
                tx_start -= min_pos
                tx_end   -= min_pos

                number_of_exons = len(es)
                pepseq = None
                if len(ee) !=  number_of_exons:
                    print "number of exon starts != number of exon ends (?) " # I'm not going there
                else:
                    # mrna
                    mrna = "".join ([seq[es[n]:ee[n]+1] for n in  range(number_of_exons)])
                    print >>outf_mrna, ">" + stable_id
                    print >>outf_mrna, mrna.lower()

                    if biotype != 'protein_coding': continue

                    # peptide
                    # to be coding, exon cannot end before the translation start, and it cannot start after the translation end
                    coding_exon_idx = [ n for  n in range(number_of_exons) if  ee[n] >= tx_start and es[n] <= tx_end]
                    dna = "".join ([seq[max(es[n],tx_start): min(ee[n],tx_end)+1] for n in  coding_exon_idx])
                    dnaseq = Seq (dna, generic_dna)
                    if strand < 0:   dnaseq = dnaseq.reverse_complement()

                    pepseq = str(dnaseq.translate())
                    print >>outf_pep, ">" + stable_id
                    print >>outf_pep, pepseq
                    if verbose:
                        print
                        print line
                        print "exon starts: ", es
                        print "exon ends: ", ee
                        print "translation start, end: ", tx_start, tx_end
                        print coding_exon_idx
                        print pepseq, "\n"
                        print


            outf_mrna.close()
            outf_pep.close()
            inf.close()

    return True


#########################################
if __name__ == '__main__':
    main()
