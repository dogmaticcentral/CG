#!/usr/bin/python -u

# to speed up things make index:  (int 01_maf_somatic_mutations_table.sql)
# create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position)
# note: I am assuming then that no tumot will have the exact same mutation in both alleles, for the simple reason that I do nto see
# how would the thwo entried in the database then be distinguished from one another
# (rather if tumor_seq_allele1 == tumor_seq_allele2 != match_norm_seq_allele1  and tumor_seq_allele2 != match_norm_seq_allele2
# then I have somehting like that)


import MySQLdb
from   tcga_utils.mysql   import  *
from math  import sqrt
import os
from  scipy import stats
import  numpy


use_normal_tissue = True

#########################################
def check_barcode(cursor, barcode):

    # is this primary tumor
    elements = barcode.split('-')
    source  = elements[3][:-1]
    # source can signa additional or metastatic tumors from the same patient
    # to keep our life simple we'll just stick to primary tumors
    # indicated by source code 01, 03, 08, or 09
    if use_normal_tissue:
        if  not source in ['10', '11', '12', '14']:
            #print "non-primary:", barcode
            return False
    else:
        if  not source in ['01', '03', '08', '09']:
            print "non-primary:", barcode
            return False

    # are there any other annotations (usually indicating something is off with the sample
    switch_to_db(cursor, 'tcga_meta')
    qry = "select * from annotation where item_barcode='%s'" % barcode
    rows = search_db(cursor, qry)
    if rows:
        print "annotation found:"
        print rows
        return False

    return True

#########################################
def bad (cursor, symbol):
    if symbol[:3] == 'LOC':  return True
    if symbol[0] == 'C' and 'ORF' in symbol: return True
    if symbol[:4]== 'KIAA': return True
    switch_to_db(cursor, 'baseline')
    qry = "select locus_type from hgnc_id_translation where approved_symbol = '%s'" % symbol
    rows = search_db(cursor, qry)
    if not rows: return  True
    if not rows[0][0] == 'gene with protein product': return True

    return False

#########################################
def  blurb (symbol, description, comment, outf):
    [nobs, [min,max], mean, variance, skewness, kurtosis, left, right]  = description
    print >> outf,"===================================="
    print >> outf,"%s " % symbol
    print >> outf,"pts: %4d    min = %6.2f   max = %10.2f    "   %  ( nobs, min, max)
    print >> outf,"mean = %10.2f   stdev = %6.2f   skew = %5.2f  kurt = %8.2f"   %  (mean, sqrt(variance), skewness, kurtosis)
    print >> outf,"fract left = %5.3f     fract right  = %5.3f " %  (left, right)
    print >> outf, comment

#########################################
def blurb_simple (symbol, description):
    [nobs, [min,max], mean, variance, skewness, kurtosis]  = description
    print "%15s " % symbol,
    print"%4d  %6.2f  %10.2f   " % (nobs, min, max),
    print "%10.2f   %6.2f   %5.2f   %8.2f" % (mean, sqrt(variance), skewness, kurtosis)


#########################################
def is_bimodal(histogram):
    for i in range (1, len(histogram)):

        val = histogram[i]
        left_max = max(10,val)
        left_max_found = False
        for j in range(0,i):
            if histogram[j] > left_max:
                left_max = histogram[j]
                left_max_found = True
                break

        right_max = max(10,val)
        right_max_found = False
        for j in range(i+1, len(histogram)):
            if histogram[j] > left_max:
                right_max = histogram[j]
                right_max_found = True
                break

        if left_max_found  and  right_max_found:
            return True

    return False
#########################################
def find_ref_values(cursor, values):

    ref_values  = {}   # for the genes showing a nice normal distribution
    for symbol, val_array in values.iteritems():
        if bad (cursor, symbol): continue
        if not val_array: continue

        description = stats.describe(val_array)
        [nobs, [min,max], mean, variance, skewness, kurtosis]  = description

        if mean < 3: continue

        if len(val_array) > 20:
            [teststat, pval] = stats.normaltest(val_array)
        else:
            teststat  = 100


        if teststat >= 4: continue

        if mean < 10: continue

        ref_values[symbol] = description

        continue # below is some descriptive output

        in_left_tail  = float (len([x for x in val_array if x < mean-2*stdev ]))/len(val_array)
        in_right_tail = float (len([x for x in val_array if x > mean+2*stdev ]))/len(val_array)
        description += (in_left_tail, in_right_tail)

        blurb(symbol, description, "normal candidate", sys.stdout);
        [hist_numpy, low, binsize, extrapoints] = stats.histogram (val_array)
        histogram = hist_numpy.tolist()
        i = 0
        bin_prev = low
        for val in histogram[:-1]:
            print " %5d %5d " % (int(bin_prev), int(bin_prev+binsize)),
            print " %5d " %  val
            i += 1
            bin_prev += binsize

        print

    return ref_values

#########################################
def process_data_set (cursor, parent_dirs, files_in_data_dir, outf):

    values = {}
    symbol_ct = 0
    for parent_dir in parent_dirs:
        # the magical formula in check_barcode gives me the barcode from the name of the file
        filtered_files = [x for x in files_in_data_dir[parent_dir] if check_barcode(cursor, x.split('.')[1]) ]
        print "parent dir:", parent_dir, "number of files:", len(filtered_files)
        for data_file in filtered_files:

            full_path = "/".join([parent_dir, data_file])
            inff = open(full_path)
            for line in inff.readlines( )[1:]:
                if line[0] == '?': continue
                line = line.rstrip()
                fields = line.split ('\t')
                if len(fields) != 4: continue # I don't know what this is in that case
                fields_clean = [x.replace("'", '') for x in fields]
                symbol = fields_clean[0].split('|')[0].upper()
                if not values.has_key(symbol):
                    values[symbol] = []
                rpkm = float (fields_clean [-1])
                #if rpkm > 2000: continue # what the hell was this?
                values[symbol].append(rpkm)


    # pass one - find near normals
    ref_values = find_ref_values(cursor, values)

    #for symbol, description in ref_values.iteritems():
    #    blurb_simple(symbol, description)

    return ref_values

    #print "normal:", normal
    #mv = [x for x in mean_values if x < 1000 ]
    #rint "mean values tot:", len(mean_values), " gt 1000: ",  len( [x for x in mean_values if x >= 1000 ])
    #[hist_numpy, low, binsize, extrapoints] = stats.histogram (mv)
    #histogram = hist_numpy.tolist()
    #if True:
    #   i = 0
    #    bin_prev = low
    #    for val in histogram[:-1]:
    #        print " %5d %5d " % (int(bin_prev), int(bin_prev+binsize)),
    #        print " %5d " %  val
    #        i += 1
    #        bin_prev += binsize
    #
    #    print


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["BRCA", "COAD","GBM","KIRC","KIRP","LAML","LGG","LUAD","LUSC","OV","REA","UCEC"]
    #db_names  = ["BRCA", "KIRC"]

    switch_to_db (cursor, 'tcga_meta')

    outf = sys.stdout
    ref_values= {}
    for db_name in db_names:

        print " ** ", db_name

        table = 'gene_expression'

        db_dir  = '/Users/ivana/databases/TCGA/%s/Expression_Genes' % db_name

        parent_dirs = []
        files_in_data_dir = {}
        for parent_dir, subdir_list, file_list in os.walk(db_dir):
            data_files = filter (lambda x: x[-3:] == 'txt', file_list)
            if data_files:
                parent_dirs.append(parent_dir)
                files_in_data_dir[parent_dir] = data_files
        if not parent_dirs:
            print "no data found for db_name (?)"
        else:
            ref_values[db_name] = process_data_set (cursor, parent_dirs, files_in_data_dir, outf)

    common_reference_genes = ([])
    for db_name in db_names:
        if not db_name in ref_values.keys() or not len(ref_values[db_name]): continue
        gene_set = set(ref_values[db_name].keys())
        if len(common_reference_genes) == 0 :
            common_reference_genes = gene_set
        else:
            common_reference_genes &= gene_set
            print
            print "after ", db_name, "common genes: ", len (common_reference_genes)

    for symbol in common_reference_genes:
        print
        for db_name in db_names:
            print " %10s" % db_name,
            blurb_simple (symbol, ref_values[db_name][symbol])

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
