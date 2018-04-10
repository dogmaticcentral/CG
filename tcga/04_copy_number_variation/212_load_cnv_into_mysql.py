#!/usr/bin/python -u

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

import os
from os import walk
import MySQLdb
from tcga_utils.mysql   import  *
from tcga_utils.ensembl  import  *
from time      import  time

import collections

#########################################
def read_cancer_names ():
    full_name= {}
    inf = open ("/Users/ivana/pypeworks/tcga/cancer_names.txt", "r")
    for line in inf:
        line   = line.rstrip() 
        field = line.split ("\t")
        if field[0] == 'READ':
            field[0] = 'REA'
        full_name[field[0]] = field[1]
    inf.close()

    return full_name


#########################################
def store (cursor, stable_id, source, seg_mean):

    # there must be no mistake, bcs we are not checking for duplicates (too slow)
    qry = "insert into cnv_snp (source_set, gene_stable_id, nocnv_seg_mean) values "
    qry += "( '%s',  '%s',  %f)" % (source, stable_id, seg_mean)
    rows   = search_db (cursor, qry)
    if (rows):
        rows   = search_db (cursor, qry, verbose=True)
        exit(1)
        return False
   

    return True

#########################################
def store_genes_seen(cursor_tcga, gene_stable, genes_seen, samples_per_gene, fold_change_per_gene, num_probes_per_gene, range_per_gene):

    for gene_id in genes_seen:
        samples = samples_per_gene[gene_id].rsplit(';')
        if len(samples) <= 1: continue

        fold_changes =  fold_change_per_gene[gene_id].rsplit(';')
        num_pbs      =  num_probes_per_gene[gene_id].rsplit(';')
        ranges       =  range_per_gene[gene_id].rsplit(';')
        for i in range (len(samples)):
            
            start, end = ranges[i].rsplit("-")
            start = int(start)
            end   = int(end)


        fixed_fields  = {'gene_stable_id': gene_stable[gene_id]}
        update_fields = {'source_ids': samples_per_gene[gene_id], 
                         'log_fold_changes': fold_change_per_gene[gene_id], 
                         'num_probes': num_probes_per_gene[gene_id], 
                         'ranges': range_per_gene[gene_id]}
        store_or_update (cursor_tcga, 'cnv_snp', fixed_fields, update_fields, verbose=False)
  
    return True

#########################################
def shutdown (cursor, cursor2, db, db2):
    cursor2.close()
    cursor.close()
    db.commit()
    db2.commit()
    db2.close()
    db.close()

#########################################
def is_sane(sample_file, seqregion_name2id, genes_per_region, gene_coordinates):

    inf = open(sample_file)
    header = []
    chrom_idx = 1
    from_idx  = 2
    to_idx    = 3
    within = 0

    for line in inf:
        if not header:
            header = line.rstrip().replace (" ", "").split("\t")
            for must_have in ['Chromosome', 'Start', 'End', 'Segment_Mean']:
                if not must_have in header:
                    print must_have, "field not found in header of ", sample_file
                    exit(1)
            chrom_idx = header.index('Chromosome')
            from_idx  = header.index('Start')
            to_idx    = header.index('End')
            seg_idx   = header.index('Segment_Mean')
            continue

        fields = line.rstrip().replace (" ", "").split("\t")
        [chromosome, start, end, seg_mean] =  [fields[chrom_idx], fields[from_idx], fields[to_idx], fields[seg_idx]]
        # picking the cutoff:
        # pow (2, 0.1) = 1.07,  pow (2, -0.1) = 0.93
        # pow (2, 0.2) = 1.15,  pow (2, -0.1) = 0.87
        if ( abs(float(seg_mean)) < 0.1) : continue
        seq_region_id = seqregion_name2id[chromosome]
        for gene_id in genes_per_region[seq_region_id]:
            [seq_region_start, seq_region_end] = gene_coordinates[gene_id]
            if ( seq_region_end < int(start) or  int(end) < seq_region_start): continue
            within += 1
                    
    inf.close()

    return (within<5000)


#########################################
def main():
    
    tcga_dir = '/Users/ivana/databases/TCGA'

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS", "UVM"];
    full_name = read_cancer_names ()

    db     = connect_to_mysql()
    cursor      = db.cursor()
    db2     = connect_to_mysql()
    cursor_tcga = db2.cursor()

    db_name = 'homo_sapiens_core_75_37'
    switch_to_db (cursor, db_name)

    start_time = time()
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1, ref_only=True)
    print "genes to consider: ", len(gene_ids), "   time: %8.3f" % (time()-start_time);

    start_time = time()
    genes_per_region = {}
    gene_coordinates = {}
    gene_stable = {}
    for gene_id in gene_ids:
        qry = "select seq_region_id, seq_region_start, seq_region_end from gene where gene_id = %d " % gene_id
        rows = search_db (cursor, qry)
        if len(rows) > 1: 
            print "bleep 1 ?  seq_region_id"
            exit(1)
        if len(rows[0]) != 3: 
            print "bleep 2 ?  seq_region_id"
            exit(1)
        [seq_region_id, seq_region_start, seq_region_end] = rows[0]
        if not genes_per_region.has_key (seq_region_id): genes_per_region[seq_region_id] = []
        genes_per_region[seq_region_id].append(gene_id)
        gene_coordinates[gene_id] = [seq_region_start, seq_region_end] 
        gene_stable[gene_id] = gene2stable (cursor, gene_id)
        
        
    print "all coords found  and grouped,  time: %8.3f" % (time()-start_time);

    seqregion_name2id = get_seq_region_ids(cursor)

    for db_name in db_names:
        print "######################################"
        print db_name, full_name[db_name]
        start_db = time()
        switch_to_db (cursor_tcga, db_name)

        path = tcga_dir + "/" + db_name + "/" +  "CNV_SNP_Array"

        if not os.path.exists(path):
            print path, "not found"
            continue

        diferent_samples_count  = 0
        uniq_samples = []
        sample_files = []
        for (dirpath, dirnames, filenames) in walk(path):
            for dirname in dirnames:
                if not 'broad' in dirname: continue
                path2 = path + "/" + dirname
                for (dirpath2, dirnames2, filenames2) in walk(path2):
                    for fnm in  filenames2:
                        if not 'nocnv' in fnm or not 'hg19' in fnm: continue
                        fields = fnm.split(".")
                        sample = fields[0]
                        if sample in uniq_samples:
                            print "duplicate sample: ", path2, sample
                        else:
                            uniq_samples.append(sample)
                            sample_files.append(path2 + "/" + fnm)
                            diferent_samples_count += 1

        ct = 0
        
        genes_seen = []
        ct2source_name = {}
        samples_per_gene = {}
        fold_change_per_gene = {}
        num_probes_per_gene = {}
        range_per_gene = {}
        
        for sample_file in sample_files:

            # do sanity check for the whole file
            # if more than a, 200 genes are reported as having significant cnv, drop the whole sample
            # the 200 cutoff comes from inspection of the output of 210_*.py - more than 1/2 sample cases in that bracket
            if not is_sane(sample_file, seqregion_name2id, genes_per_region, gene_coordinates): continue

            ct += 1
            start_time = time()
            source = sample_file.split ('/')[-1]
            ct2source_name[str(ct)] = source
            
            fixed_fields  = {'source_id':ct, 'source_name':source}
            update_fields = {}
            store_or_update (cursor_tcga, 'cnv_source_ids', fixed_fields, update_fields, verbose=False)
            if not ct%100: 
                print source, ct, "out of ", len(sample_files)


            inf = open(sample_file)
            header = []
            chrom_idx = 1
            from_idx  = 2
            to_idx    = 3

            for line in inf:
                if not header:
                    header = line.rstrip().replace (" ", "").split("\t")
                    for must_have in ['Chromosome', 'Start', 'End', 'Segment_Mean', 'Num_Probes']:
                        if not must_have in header:
                            print must_have, "field not found in header of ", sample_file
                            exit(1)
                    chrom_idx = header.index('Chromosome')
                    from_idx  = header.index('Start')
                    to_idx    = header.index('End')
                    seg_idx   = header.index('Segment_Mean')
                    num_idx   = header.index('Num_Probes')
                    continue

                fields = line.rstrip().replace (" ", "").split("\t")
                
                [chromosome, start, end,  num_probes, seg_mean] =  [fields[chrom_idx], fields[from_idx], fields[to_idx], fields[num_idx], fields[seg_idx]]
                # picking the cutoff:
                # pow (2, 0.1) = 1.07,  pow (2, -0.1) = 0.93
                # pow (2, 0.2) = 1.15,  pow (2, -0.1) = 0.87
                if ( abs(float(seg_mean)) < 0.1): continue
                
                seg_mean = "%.2f" % float(seg_mean)

                if not seqregion_name2id.has_key(chromosome):
                    print "seq region not found", chromosome
                    exit(1)

                #####################
                seq_region_id = seqregion_name2id[chromosome]

                for gene_id in genes_per_region[seq_region_id]:
                    
                    [seq_region_start, seq_region_end] = gene_coordinates[gene_id]
                    if ( seq_region_end < int(start) or  int(end) < seq_region_start): continue
                    if not gene_id in genes_seen:  
                        genes_seen.append(gene_id)
                        # we are using ct as source identifier
                        samples_per_gene[gene_id]     = "%s" % ct
                        fold_change_per_gene[gene_id] = seg_mean
                        num_probes_per_gene[gene_id]  = num_probes
                        range_per_gene[gene_id]       = start + "-" + end
                    else:
                        samples_per_gene[gene_id]     += ";%s" % ct
                        fold_change_per_gene[gene_id] += ";" + seg_mean
                        num_probes_per_gene[gene_id]  += ";" + num_probes
                        range_per_gene[gene_id]       += ";" + start + "-" + end
              
            #print "\t\t time: %8.3f" % (time()-start_time)

        print db_name, "number of genes seen:", len(genes_seen), "in  %8.3f s" % (time()-start_db)
        store_genes_seen(cursor_tcga, gene_stable, genes_seen, samples_per_gene, fold_change_per_gene, num_probes_per_gene, range_per_gene)


        db.commit()
        db2.commit()

    shutdown (cursor, cursor_tcga, db, db2)



#########################################
if __name__ == '__main__':
    main()
