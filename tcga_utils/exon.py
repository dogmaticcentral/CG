
import os

#########################################################
class Exon:



    #################################################
    # the ensembl row is assumed to have been 
    # obtained by select * from exon
    def load_from_ensembl_exon(self, gene_start, gene_end, ensembl_row):
 
        self.exon_id          = ensembl_row[0]
        strand =  ensembl_row[4]
        start_on_seq_region  = ensembl_row[2]
        end_on_seq_region    = ensembl_row[3]

        if ( strand > 0 ):
            self.start_in_gene = start_on_seq_region-gene_start
            self.end_in_gene   = end_on_seq_region-gene_start
        else:
            self.start_in_gene = gene_end - end_on_seq_region
            self.end_in_gene   = gene_end - start_on_seq_region

        self.strand            = strand
        self.phase             = ensembl_row[5]
        self.end_phase         = ensembl_row[6]
        self.is_constitutive   = ensembl_row[8]
        #  known (that's the source indicator - we 
        # got it from table called 'exon')
        self.is_known          = 1
        
        
        return True

    #################################################
    # the ensembl row is assumed to have been 
    # obtained by select * from exon_pprediction
    def load_from_ensembl_prediction(self, gene_start, gene_end, ensembl_row):

        self.exon_id           = ensembl_row[0]
        strand = ensembl_row[6]
        start_on_seq_region = ensembl_row[4]
        end_on_seq_region   = ensembl_row[5]

        if ( strand > 0 ):
            self.start_in_gene = start_on_seq_region-gene_start
            self.end_in_gene   = end_on_seq_region  -gene_start
        else:
            self.start_in_gene = gene_end - end_on_seq_region
            self.end_in_gene   = gene_end - start_on_seq_region
        self.strand            = strand
        self.phase             = ensembl_row[7]
        # not known (that's the source indicator - we 
        # got it from table called 'prediction_exon')
        self.is_known          = 0
        self.is_canonical      = 0
        self.is_constitutive   = 0
        
        return True


    #################################################
    # load in from gene2exon table (select * from gene2exon)
    def load_from_gene2exon (self, gene2exon_row):

        if ( len(gene2exon_row) < 17):
            print "error loading exon: the in list must be",
            print " at least 17 elements long"
            return False
        
        self.gene_id             = gene2exon_row [1]
        self.exon_id             = gene2exon_row [2]
        self.start_in_gene       = gene2exon_row [3]
        self.end_in_gene         = gene2exon_row [4]
        self.canon_transl_start  = gene2exon_row [5]
        self.canon_transl_end    = gene2exon_row [6]
        self.exon_seq_id         = gene2exon_row [7]
        self.strand              = gene2exon_row [8]
        self.phase               = gene2exon_row [9]
        self.is_known            = gene2exon_row[10]
        self.is_coding           = gene2exon_row[11]
        self.is_canonical        = gene2exon_row[12]
        self.is_constitutive     = gene2exon_row[13]
        self.covering_exon       = gene2exon_row[14]
        self.covering_exon_known = gene2exon_row[15]
        self.analysis_id         = gene2exon_row[16]

        return True

    #################################################
    # load in from gene2exon table (select * from gene2exon)
    def load_from_novel_exon (self, table_row, table):

        if ( len(table_row) < 11):
            print "error loading exon: the in list must be",
            print " at least 10 elements long"
            return False
        
        self.exon_id               = table_row[0]
        self.gene_id               = table_row[1]
        self.start_in_gene         = table_row[2]
        self.end_in_gene           = table_row[3]
        self.maps_to_human_exon_id = table_row[4]
        self.exon_seq_id           = table_row[5]
        self.template_exon_seq_id  = table_row[6]
        self.template_species      = table_row[7]

        self.strand              = table_row [8]
        self.phase               = table_row [9]
        self.end_phase           = table_row[10]
        self.has_NNN             = table_row[11]
        self.has_stop            = table_row[12]
        self.has_3p_ss           = table_row[13]
        self.has_5p_ss           = table_row[14]
        if table=='sw_exon':
            self.is_known     =  2 
            self.analysis_id  = -1 
        else:
            self.is_known     =  3
            self.analysis_id  = -2
           
        self.is_coding           = 1
        self.is_canonical        = 0
        self.is_constitutive     = 0
        return True

    #################################################
    # print
    def __str__ (self):

        printstr = ""

        for attr, value in self.__dict__.iteritems():
            if ( not value is None):
                printstr += " %-20s    %s" % (attr,  str(value))
            else:
                printstr += " %-20s    None"  % attr
            printstr += "\n"

        return printstr

    ###################################
    # when something is defined as an exon ....
    def __init__ (self):
        
        self.exon_id               = None
        self.exon_seq_id           = None
        self.stable_id             = None
        self.gene_id               = None
        self.start_in_gene         = None
        self.end_in_gene           = None
        self.exon_seq_id           = None
        self.strand                = None
        self.phase                 = None
        self.end_phase             = None
        self.pepseq_transl_start   = None
        self.pepseq_transl_end     = None
        self.canon_transl_start    = None
        self.canon_transl_end      = None
        self.is_known              = None
        self.is_coding             = None
        self.is_canonical          = None
        self.covering_exon         = None
        self.covering_exon_known   = None
        self.maps_to_human_exon_id = None
        self.analysis_id           = None
        self.pepseq                = None
        self.template_species      = None
        self.template_exon_seq_id  = None
        self.has_NNN               = None
        self.has_stop              = None
        self.has_3p_ss             = None
        self.has_5p_ss             = None

