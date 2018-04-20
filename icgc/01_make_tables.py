#!/usr/bin/python

import MySQLdb
from icgc_utils.mysql   import  *

# icgc_donor_id	 icgc_specimen_id	icgc_sample_id
# chromosome	 chromosome_start	chromosome_end	chromosome_strand	assembly_version
# mutation_type	 reference_genome_allele	mutated_from_allele	 mutated_to_allele
# consequence_type	aa_mutation	 cds_mutation
# gene_affected	transcript_affected
#########################################
def make_mutations_table(cursor, db_name, mutations_table):

    switch_to_db (cursor, db_name)

    qry = ""
    qry += "  CREATE TABLE  %s (" % mutations_table
    qry += "     id INT NOT NULL AUTO_INCREMENT, "
    qry += "  	 hugo_symbol VARCHAR (50) NOT NULL, "
    qry += "     entrez_gene_id INT, "
    qry += "	 aa_change VARCHAR (100), "
    qry += "	 cdna_change BLOB, "
    qry += "	 chromosome VARCHAR (20) NOT NULL, "
    qry += "	 start_position INT  NOT NULL, "
    qry += "	 end_position INT NOT NULL, "
    qry += "	 strand VARCHAR (5) NOT NULL, "
    qry += "	 variant_classification VARCHAR (50) NOT NULL, "
    qry += "	 variant_type VARCHAR (20) NOT NULL, "
    qry += "	 reference_allele BLOB NOT NULL, "
    qry += "	 tumor_seq_allele1 BLOB NOT NULL, "
    qry += "	 tumor_seq_allele2 BLOB NOT NULL, "
    qry += "	 tumor_sample_barcode VARCHAR (50) NOT NULL, "
    qry += "	 sample_barcode_short VARCHAR (20) NOT NULL, "
    qry += "	 matched_norm_sample_barcode VARCHAR (50) NOT NULL, "
    qry += "	 match_norm_seq_allele1 BLOB, "
    qry += "	 match_norm_seq_allele2 BLOB, "
    qry += "	 tumor_validation_allele1 BLOB, "
    qry += "	 tumor_validation_allele2 BLOB, "
    qry += "	 match_norm_validation_allele1 BLOB, "
    qry += "	 match_norm_validation_allele2 BLOB, "
    qry += "	 verification_status VARCHAR (20), "
    qry += "	 validation_status VARCHAR (20) NOT NULL, "
    qry += "	 mutation_status VARCHAR (50) NOT NULL, "
    qry += "	 conflict BLOB, "
    qry += "	 PRIMARY KEY (id) "
    qry += ") ENGINE=MyISAM"

    rows = search_db(cursor, qry)
    print qry
    print rows

#########################################
#########################################
def main():

    db     = connect_to_mysql()
    cursor = db.cursor()

