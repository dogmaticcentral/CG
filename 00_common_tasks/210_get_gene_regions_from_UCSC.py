#!/usr/bin/python
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
from tcga_utils.mysql import *
import os
#########################################
def construct_qry (chrom, fields):
    qry = "select "
    first = True
    for field in fields:
        if first:
            first = False
        else:
            qry += ", "
        qry += "g." + field
    qry += " "
    qry += "from refGene as g, gbCdnaInfo as i, refSeqStatus as s "
    qry += "where g.chrom='%s' " % chrom
    qry += "and i.acc=g.name  and i.type='mRNA' and i.mol='mRNA' "
    qry += "and s.mrnaAcc=g.name and s.status in ('Validated', 'Reviewed')"

    return qry

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
        fields_to_dwld = ["name", "name2", "strand","txStart", "txEnd", "exonStarts",  "exonEnds"]
        target_dir = "/mnt/databases/ucsc/gene_regions/" + assembly
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        for chrom in chromosomes:
            print "downloading data for", chrom
            outf = open(target_dir + "/" + chrom+".csv", "w")
            print  >>outf,  "\t".join( fields_to_dwld )
            qry = construct_qry(chrom, fields_to_dwld)
            rows = search_db(cursor,qry)

            for row in rows:
                print  >>outf,  "\t".join( [ str(r) for r in row] )
            outf.close()

    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


