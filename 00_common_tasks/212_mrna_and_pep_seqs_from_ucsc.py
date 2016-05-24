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
        target_dir = {}
        for seq_type in ['mrna', 'pep']:
            target_dir[seq_type] = "/mnt/databases/ucsc/sequences/" + assembly + "/" + seq_type
            if not os.path.isdir(target_dir[seq_type]):
                os.makedirs(target_dir[seq_type])
        for chrom in chromosomes:
            print "downloading mrna for", chrom
            outf = open(target_dir['mrna'] + "/" + chrom+".fasta", "w")
            qry  = "select m.name, m.seq from knownGeneMrna as m, knownCanonical as g "
            qry += " where m.name = g.transcript and g.chrom = '%s'"  % chrom
            rows = search_db(cursor, qry)
            if not rows:
                print "no return for ", qry
            for row in rows:
                [name, seq] = row
                print >>outf,  ">" + name
                print >>outf, seq
            outf.close()

        for chrom in chromosomes:
            print "downloading peptides for", chrom
            outf = open(target_dir['pep'] + "/" + chrom+".fasta", "w")
            qry  = "select m.name, m.seq from knownGenePep as m, knownCanonical as g "
            qry += " where m.name = g.protein and g.chrom = '%s'"  % chrom
            rows = search_db(cursor, qry)
            if not rows:
                print "no return for ", qry
            for row in rows:
                [name, seq] = row
                print >>outf, ">" + name
                print >>outf, seq
            outf.close()

    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


