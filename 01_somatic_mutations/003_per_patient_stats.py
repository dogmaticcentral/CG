#!/usr/bin/python -u
# the fields in the tumor sample barcode are
# project - tissue source site (TSS)  - participant -
# source.vial - portion.analyte  - plate - (sequencing or characterization center)
#
#
# sample codes:
# 01	Primary solid Tumor	TP
# 02	Recurrent Solid Tumor	TR
# 03	Primary Blood Derived Cancer - Peripheral Blood	TB
# 04	Recurrent Blood Derived Cancer - Bone Marrow	TRBM
# 05	Additional - New Primary	TAP
# 06	Metastatic	TM
# 07	Additional Metastatic	TAM
# 08	Human Tumor Original Cells	THOC
# 09	Primary Blood Derived Cancer - Bone Marrow	TBM
# 10	Blood Derived Normal	NB
# 11	Solid Tissue Normal	NT
# 12	Buccal Cell Normal	NBC
# 13	EBV Immortalized Normal	NEBV
# 14	Bone Marrow Normal	NBM
# 20	Control Analyte	CELLC
# 40	Recurrent Blood Derived Cancer - Peripheral Blood	TRB
# 50	Cell Lines	CELL
# 60	Primary Xenograft Tissue	XP
# 61	Cell Line Derived Xenograft Tissue	XCL


import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *

verbose = False

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()


    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    table = 'somatic_mutations'
    total_number_of_somatic_mutations_in_TCGA = 0
    
    for db_name in db_names:

        switch_to_db (cursor, db_name)
        print 

        ############################
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        print db_name, table, ", number of entries (somatic mutations):", rows[0][0]
        total_number_of_somatic_mutations_in_TCGA += int(rows[0][0])
        ############################

        ############################
        if verbose: # so far I haven't seen an anetry without tumor sample barcode, but you never know with TCGA
            qry  = "select count(1) from somatic_mutations "
            qry += "where not tumor_sample_barcode like 'TCGA%'"
            rows = search_db(cursor, qry)
            print "number of entries without tumor barcode:", rows[0][0]
        ############################

        ############################
        if verbose: print "sorting per-sample data ..."
        qry  = "select distinct  sample_barcode_short from somatic_mutations "
        rows = search_db(cursor, qry)
        print "number of different samples:", len(rows)
        # now the number of samples split into different vials
        mult_vials_or_info = ""
        mvcount = 0
        sources_per_tumor = {}
        multiple_sources = False
        for  row in rows:
            sample_barcode_short = row[0]
            qry = "select distinct  tumor_sample_barcode from somatic_mutations "
            qry += "where  sample_barcode_short = '%s'"  % sample_barcode_short
            rows2 = search_db(cursor, qry)
            if len(rows2)>1:
                mvcount += 1
                mult_vials_or_info +=  "\t" + str(sample_barcode_short) +  ": " + str(rows2) + "\n"
            # sample_barcode_short is tissue source site (TSS)  - participant - source.vial
            fields = sample_barcode_short.split('-')
            tumor_id =   '-'.join(fields[:2])
            source  = fields[2][:2]
            if not sources_per_tumor.has_key(tumor_id):
                sources_per_tumor[tumor_id] = []
            else:
                multiple_sources = True
            sources_per_tumor[tumor_id].append(source)

        if mvcount:
            print "samples from multiple vials, portion, plate, or different seq or characterization center (%d)" % mvcount
            print "this should not really be a problem - we are removing duplicate entries (same sample, same position on the genome)"
            if verbose: print mult_vials_or_info
        if multiple_sources:
            print "patients contributing multiple sources - potentially interesting, but how do we handle this?"
            # (for now we'll solve it by going back to filling the database, and dropping all samples
            # that are not 01, 03, 08 or 09 (see table on top of the file)
            for patient, sources in sources_per_tumor.iteritems():
                if len(sources) < 2: continue
                print "\t", patient, sources
    print
    print 'total number of somatic mutations in TCGA:', total_number_of_somatic_mutations_in_TCGA
    print
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

