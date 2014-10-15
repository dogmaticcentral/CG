#!/usr/bin/python -u

import MySQLdb
from sets import Set
from   tcga_utils.mysql   import  *
import commands, re

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["COAD",  # after this we go alphabetically
                 "ACC", "BLCA", "BRCA", "CESC",  "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]


    position_pattern   = re.compile ('[pc]\.\D*(\d+)\D*')
    from_to_pattern    = re.compile ('[pc]\.[\d\>]*(\D*)[\d\>]+(\D*)')

    for db_name in db_names:

        print " ** ", db_name
        switch_to_db (cursor, db_name)

        table = 'somatic_mutations'

        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name
            exit(1)

        qry = "select hugo_symbol, aa_change, cdna_change from somatic_mutations "
        rows = search_db (cursor, qry)

        if not rows:
            print "no somatic mutations found"
            continue
        does_not_compute = 0
        computes = 0
        for row in rows:
            [hugo_symbol, aa_change, cdna_change] = row
 
            aa_position_match   = position_pattern.match (aa_change)
            cdna_position_match = position_pattern.match (cdna_change)

            aa_from_to_match   = from_to_pattern.match (aa_change)
            cdna_from_to_match = from_to_pattern.match (cdna_change)

            if not aa_position_match:
                #print "\t\t no aa position match"
                pass
            elif not cdna_position_match:
                #print "\t\t no cdna position match"
                pass
            elif not aa_from_to_match:
                #print "\t\t no aa from to  match"
                pass
            elif not cdna_from_to_match:
                #print "\t\t no cdna from to  match"
                pass
            elif not  (len(aa_from_to_match.group(1))==1  and len(aa_from_to_match.group(2)) == 1) :
                #print "\t ", hugo_symbol, aa_change, cdna_change
                pass # let's look for point mutations for now
            elif aa_from_to_match.group(2) == '*':
                pass
            elif not  (len(cdna_from_to_match.group(1))==1  and len(cdna_from_to_match.group(2)) == 1) :
                # double mutations seem to be quite common in SKCM, but even there there are only
                # 982 compred to 263862 SNPs (thus, we ignore them for now)
                pass

            else:
                
                position_protein =  int (aa_position_match.group(1))
                position_dna     =  int (cdna_position_match.group(1))
                # are the protein and dna positions related as expected
                pos_computes = (position_dna%3)  and (position_protein == position_dna/3+1) \
                               or  position_protein == position_dna/3;
                if not pos_computes:
                    print "\t ", hugo_symbol, aa_change, cdna_change, 
                    does_not_compute += 1
                else:
                    computes +=1
            


                #print "\t\t ", aa_from_to_match.group(1), cdna_from_to_match.group(1)
                #print "\t\t ", position_protein, position_dna, position_dna/3, position_dna%3
                #print "\t\t ", aa_from_to_match.group(2), cdna_from_to_match.group(2)
                pass
        print "computes: ", computes, "does not compute: ", does_not_compute



#########################################
if __name__ == '__main__':
    main()
