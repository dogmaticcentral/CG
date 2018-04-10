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

# I am not sure how to check that the level3 data for CNV make sense
# the amount of info is really sparse. The columns:
# Sample	Chromosome	Start	End	Num_Probes	Segment_Mean
# and no cross-refernces are provided

# I must admit I cringe at the prospect of going to level2, and I am
# not even sure that it is available

# one single thing that can be done is running through all data and making sure that
# the coordinates exist on the given chromosome


# as it turns out,  the coordinates have changed little enough so that both hg18 and hg38 give me 
# about the same output

import os
from os import walk
import MySQLdb
from   tcga_utils.mysql   import  *

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
def  get_seq_region_ids(cursor):

    seqregion_name2id = {}

    qry  = "select seq_region.name, seq_region.seq_region_id from seq_region, "
    qry += "coord_system where coord_system.coord_system_id = seq_region.coord_system_id "
    qry += "and coord_system.attrib='default_version' "
    rows = search_db (cursor, qry)
    if ( not rows):
        print "bleep?"
        exit(1)
    for row in rows:
        seqregion_name2id[ str(row[0]) ] = str(row[1])
    return seqregion_name2id

#########################################
def main():
    
    tcga_dir = '/Users/ivana/databases/TCGA'

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS", "UVM"];
    full_name = read_cancer_names ()

    db     = connect_to_mysql()
    cursor = db.cursor()

    switch_to_db (cursor, 'homo_sapiens_core_75_37')
    #switch_to_db (cursor, 'homo_sapiens_core_78_38')
    seqregion_name2id = get_seq_region_ids(cursor)

    for db_name in db_names:
        print "######################################"
        print db_name, full_name[db_name]
        
        path = tcga_dir + "/" + db_name + "/" +  "CNV_SNP_Array"

        if not os.path.exists(path):
            print path, "not found"
            continue

        print path, "ok"
        
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

        print "number of different samples ", diferent_samples_count
        for sample_file in sample_files:
            inf = open(sample_file)
            header = []
            chrom_idx = 1
            from_idx  = 2
            to_idx    = 3
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
                    continue
                fields = line.rstrip().replace (" ", "").split("\t")
                [chromosome, start, end] =  [fields[chrom_idx], fields[from_idx], fields[to_idx]]
                if not seqregion_name2id.has_key(chromosome):
                    print "seq region not found", chromosome
                    exit(1)
                else:
                    pass
                    #print seqregion_name2id[chromosome], chromosome, start, end
                # now  check that the region exists on this chromosome as advertized
                seq_region_id = seqregion_name2id[chromosome]
                #qry  = "select gene_id, seq_region_id, seq_region_start, seq_region_end "
                #qry += "from gene  where  seq_region_end >= %d " % int(start)
                qry  = "select gene_id, seq_region_id, seq_region_start, seq_region_end, description "
                qry += "from gene where  biotype='protein_coding' and seq_region_id = %d" %  int(seq_region_id)
                qry += " and seq_region_end >= %d " % int(start)
                qry += " and  seq_region_start <= %d " % int(end)
                rows = search_db (cursor, qry)
                if ( not rows):
                    print "no corresponding region"
                else:
                    print '-'*23
                    print line
                    for row in rows:
                        print row, start, end
                    exit(1)
                    pass
            exit(1)
        #exit(1)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
