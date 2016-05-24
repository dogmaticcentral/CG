#!/usr/bin/python
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
from tcga_utils.utils import make_named_fields
from tcga_utils.ucsc import *

import os

#########################################
def make_canonical_transcript_table (cursor, db_name, tbl):
    switch_to_db (cursor, db_name)

    qry = ""
    qry += "  CREATE TABLE  %s (" % tbl
    qry += "     id INT NOT NULL AUTO_INCREMENT, "
    qry += "  	 transcript_id VARCHAR (50), "
    qry += "  	 protein_id VARCHAR (50), "
    qry += "  	 hugo_symbol VARCHAR (500), "
    qry += "	 strand VARCHAR (5) NOT NULL, "
    qry += "	 tx_start INT  NOT NULL, "
    qry += "	 tx_end INT NOT NULL, "
    qry += "	 exon_starts BLOB, "
    qry += "	 exon_ends BLOB, "
    qry += "	 mrna MEDIUMBLOB, "
    qry += "	 protein BLOB, "
    qry += "	 PRIMARY KEY (id) "
    qry += ") ENGINE=MyISAM"
    rows = search_db(cursor, qry)
    print qry
    print rows

    qry = "create index range_idx on %s (tx_start, tx_end)" % tbl
    rows = search_db(cursor, qry)
    print qry
    print rows

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
         for chrom in chromosomes:
            table_name = "canonical_transcripts_%s_%s" % (assembly, chrom)
            if ( check_table_exists (cursor, db_name, table_name)):
                print table_name, " found in ", db_name
                #qry = "drop table %s" % table_name
                #search_db(cursor, qry)
            else:
                make_canonical_transcript_table (cursor, db_name, table_name)
    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


