#!/usr/bin/python -u
# store the meta info about the maf files: name and the reference genome,  for now

import os.path
from tcga_utils.mysql   import  *
from tcga_utils.utils   import  *
import commands


#########################################
def  find_reference_genome (maffile):

    ref_gen = ""
    # find 10 examples from each chromosome
    inff = open(maffile, "r")
    line_ct = 0
    for line in inff:
        line_ct += 1
        if  line_ct%1000: continue
    inff.close()
    # try matching using ucsf server
    # use the assembly with the smallest amount of failures
    # if both assemblies fail in more htan 10% of cases - abort
    return ref_gen


#########################################
def store_meta_info (cursor,  maf_name, ref_gen):

    return

#########################################
def check_headers (maffile, required_fields, expected_fields):

    # required_fields are the absolute minimum we need to 
    # reconstruct the mutation - that they are missing  should not happen at all    
    inff = open(maffile, "r")
    missing_fields = []
    headerline = ""
    for line in inff:
        if line.isspace(): continue
        if line[0]=='#': continue
        headerline = line.rstrip()
        break # the first line that is not the comment should be the header
    inff.close()
    
    header_fields  = process_header_line(headerline)
    if len(header_fields) == 0:
        return ["fail", "no header found"]
    missing_fields = filter (lambda x: not x in header_fields, required_fields)
    if len(missing_fields) > 0:
        return ["fail", "missing fields:  " + " ".join(missing_fields)]

    missing_fields = filter (lambda x: not x in header_fields, expected_fields)
    if len(missing_fields) > 0:
        return ["warn", "expected fields:  " + " ".join(missing_fields)]
   
    return ["pass", ""]



##################################################################################
# checking for the following, as seen in
# broad.mit.edu_LIHC.IlluminaGA_DNASeq_automated.Level_2.1.0.0/
# An_TCGA_LIHC_External_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.maf 
#           273933        RPL5       Frame_Shift_Del         p.K270fs 
#           273933        RPL5       Frame_Shift_Del         p.K277fs 
#           273933        RPL5       Frame_Shift_Del         p.R279fs 
#           273933        RPL5       Frame_Shift_Del         p.Q282fs 
##################################################################################
def check_health (maffile):
    
    inff = open(maffile, "r")
    missing_fields = []
    headerline = ""
    for line in inff:
        if line.isspace(): continue
        if line[0]=='#': continue
        headerline = line.rstrip()
        break # the first line that is not the comment should be the header
    
    header_fields  = process_header_line(headerline)
    if len(header_fields) == 0:
        # though we shold have discovered this previously ....
        return ["fail", "no header found"]
    
    variantclass_index = header_fields.index('variant_classification')
    startpos_index     = header_fields.index('start_position')
    sample_barcode_idx = header_fields.index('tumor_sample_barcode')
    start_posns = {}
    for line in inff:
        if not 'Frame_Shift_Del' in line: continue
        field = line.split ("\t");
        if not field[variantclass_index] ==  'Frame_Shift_Del': continue
        sample_barcode = field[sample_barcode_idx]
        if not start_posns.has_key(sample_barcode): start_posns[sample_barcode] = []
        start_posns[sample_barcode].append(int(field[startpos_index]))
    inff.close()

    number_of_samples_with_stutter_count = {}
    for sample_barcode in  start_posns.keys():
        prev = None
        count = 0
        for sp in start_posns[sample_barcode]:
            if prev and sp-prev < 5:
                count += 1
            prev = sp
        if not count in number_of_samples_with_stutter_count.keys(): number_of_samples_with_stutter_count[count] = 0
        number_of_samples_with_stutter_count[count] += 1
        
    #for stutter_count in sorted( number_of_samples_with_stutter_count.keys()):
    #    print "\t", stutter_count, number_of_samples_with_stutter_count [stutter_count]

    # I am not really sure what to do with this, so for now I will
    # only store it with a warning
    bad_counts = filter (lambda x: x>100, number_of_samples_with_stutter_count.keys())
    if len(bad_counts)>0:
        return  ["warn", "%d samples have more that 100 stutter frameshift points"%len(bad_counts)]
        
    return ["pass", ""]



##################################################################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()
    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL",  "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",  "MESO", "OV",   "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
    
    
    for db_name in db_names:
        # check db exists
        qry = "show databases like '%s'" % db_name
        rows = search_db(cursor, qry)
        if not rows:
            print db_name, "not found"
            continue
        
        print " ** ", db_name
        switch_to_db (cursor, db_name)
                
        db_dir  = '/mnt/databases/TCGA/' + db_name
        if not  os.path.isdir(db_dir):
            print "directory " + db_dir + " not found"
            exit(1)
            
        ret       = commands.getoutput('find ' + db_dir + ' -name "*.maf"')
        maf_files = ret.split('\n')
        required_fields = get_required_fields()
        expected_fields = get_expected_fields (cursor, db_name, "somatic_mutations")
    
        for maffile in maf_files:
            print '\n\t processing:', bare_filename
            ref_genome = ""
            bare_filename =  maffile.split('/')[-1]
            overall_diagnostics = []

            if "automated" in bare_file.lower():
                overall_diagnostics.append(["warn", "automated"])
            elif  "curated" in bare_file.lower():
                overall_diagnostics.append(["warn", "curated"])
            
            if os.path.getsize(maffile)==0:
                diagnostics = ["fail", "file empty"]
                ref_genome  = ""
                #store (cursor, barefile, diagnostics, ref_genome)
                print "\t storing: ",  [diagnostics], ref_genome
                continue
            
            # check if the file contains  all the info we need and hope to have
            diagnostics = check_headers (maffile, required_fields, expected_fields)
            print "\t diag: ", [diagnostics]
            if diagnostics[0] == "fail":
                ref_genome  = ""
                #store (cursor, barefile, diagnostics, ref_genome)
                continue
            elif diagnostics[0] == "warn":
                overall_diagnostics.append(diagnostics)

            
            # I am aware of one way in which the file can be corrupt, so I am checking for it
            diagnostics = check_health (maffile)
            print "\t diag: ", [diagnostics]
            if diagnostics[0] == "fail":
                ref_genome  = ""
                #store (cursor, barefile, diagnostics, ref_genome)
                continue
            elif diagnostics[0] == "warn":
                overall_diagnostics.append(diagnostics)

            ref_gen = find_reference_genome (maffile)
            if not ref_gen:
                diagnostics = ["fail", "reference genome could not be determined"]
                overall_diagnostics.append(diagnostics)
                # if there is no aa_info, and no reference genome, move on to the next file
                missing = filter (lambda x: "missing" in x[1] and "aa_change" in  x[1], overall_diagnostics)
                if len(missing) > 0:
                    continue
            else:
                ref_genome = ref_gen
                diagnostics = ["pass", ""]

            #store (cursor, barefile, diagnostics, ref_genome)
            
            
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
