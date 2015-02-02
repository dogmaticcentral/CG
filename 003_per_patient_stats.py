#!/usr/bin/python -u

import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from scipy import stats

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    
    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS", "UVM"];
    table = 'somatic_mutations'

    
    for db_name in db_names:

        switch_to_db (cursor, db_name)
        print 

        ############################
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        print db_name, table, ", number of entries:", rows[0][0]
        ############################

        ############################
        qry  = "select count(1) from somatic_mutations "
        qry += "where not tumor_sample_barcode like 'TCGA%'"
        rows = search_db(cursor, qry)
        print "number of entries without tumor barcode:", rows[0][0]
        ############################

        ############################
        print "sorting per-patient data ..."
        uniq_patients = []
        tbarcodes_per_patient = {}
        qry  = "select distinct  tumor_sample_barcode from somatic_mutations "
        rows = search_db(cursor, qry)
        for  row in rows:
            tbarcode = row[0]
            # the fields are 
            # project - tissue source site (TSS)  - participant -
            # source.vial - portion.analyte  - plate - (sequencing or charcterization center)
            fields = tbarcode.split('-')
            patient = '-'.join(fields[1:3])
            if not patient in  uniq_patients:
                uniq_patients.append(patient)
                tbarcodes_per_patient[patient] = []
            tbarcodes_per_patient[patient].append(tbarcode)
        print "number of different patients:", len(uniq_patients)
        ############################

        ############################
        count = 0
        for patient in uniq_patients:
            if len( tbarcodes_per_patient[patient] ) > 1: count += 1
        print "number of patients with more than one sample:", count

        if count > 1:
            print "number of mutations in alternative samples:"
            for patient in uniq_patients:
                if len( tbarcodes_per_patient[patient] ) == 1: continue
                print patient
                for tbarcode in tbarcodes_per_patient[patient]:
                    qry = "select count(1) from somatic_mutations where tumor_sample_barcode = '%s'" % tbarcode
                    rows = search_db (cursor, qry)
                    if not rows: continue
                    print "\t ", tbarcode, rows[0][0]
       
        continue
        ############################

        ############################
        # sort by the number of entries
        print "number of different samples per patient:"
        sample_sizes = []
        patients_sorted = sorted(uniq_patients.keys(), key= lambda x: len( uniq_patients[x]) )
        for patient in patients_sorted:
            samples = uniq_patients[patient]
            print patient, len(samples)
            uniq_samples = {}
            for sample in samples:
                if not  uniq_samples.has_key(sample):
                    uniq_samples[sample] = 0
                    uniq_samples[sample] += 1
                    if len(uniq_samples) >1:
                        pass
                        #for sample in  uniq_samples.keys():
                        #    print "\t", sample, uniq_samples[sample]
                    else:
                        sample_sizes.append(len(samples))

        print
        print "histogram of the number of mutations"
        [histo, low_range, binsize, extrapoints] = stats.histogram(sample_sizes, numbins=500)
        prev = low_range
        for entry in histo:
            if entry: print " %6d  %6d   %5d " % ( int(prev), int(prev+binsize), int(entry))
            prev += binsize

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

