#!/usr/bin/python
#
#
# downloading from http to a file (using urllib):
# http://stackoverflow.com/questions/19602931/basic-http-file-downloading-and-saving-to-disk-in-python
#

#Archive Name Format 
# from https://wiki.nci.nih.gov/display/TCGA/TCGA+Data+Archives#TCGADataArchives-NamingConventions
#<domain>_<disease study>.<platform>.<archive type>.<serial index>.<revision>.<series>
#  the serial index is a number that uniquely identifies an independent data set from a particular experiment. 
#
#There is no overlap of data files between archives of differing serial numbers. 
#A numbering is entirely up to the data submission center. 
#In general, BCRs use a serial index equivalent to a batch number while other center 
#types start serial index series from 1.
#
#A revision number can indicate the number of times an archive has been revised (starting from 0) and submitted to the DCC. However, the only requirement for revision numbers is that the revision number of the new archive is to be higher than that of the archive being replaced. Files that have been changed or added are captured in the changes and additions files, respectively
#
#  the series number should always be 0.
#
# note: run twice to get the unzipped files

import urllib, os, subprocess

#########################################
def check_and_make(path):
     if not os.path.exists(path):
         print path, "not found, making", path
         os.makedirs(path)
     elif not os.path.isdir(path):
         print ">>> PROBLEM <<<<< ", path,  "is not directory ?"
         exit(1)

#########################################
def main():
    
    root = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor'
    
    local_dir = '/Users/ivana/databases/TCGA'

    file_names = open('snp_files.txt')
    all_files = {}
    for line in [x.rstrip().replace(' ','')  for x in file_names]:
        [directory, filename] = line.split(':')
        if directory[-1] == '/':   directory = directory[:-1]
        if not all_files.has_key(directory): all_files[directory] = []
        all_files[directory].append(filename)
    file_names.close()

    newest_revision = {}
    dir_for_data = {}
    for directory, files in all_files.iteritems():
        newest_revision[directory] = {}
        for fnm in files:
            fields = fnm.split ('.')
            serialno = fields[-5]
            revision = int(fields[-4])
            exp_identifier = '.'.join(fields[:-4])
            #print "\t", exp_identifier,  revision
            if not newest_revision[directory].has_key(exp_identifier) \
               or newest_revision[directory][exp_identifier] < revision:
                newest_revision[directory][exp_identifier] = revision
            #I'm not going to do do anything about the possbility that the
            # same exp data ended up in different directories - I'm just going to crash in that case
            if not dir_for_data.has_key(exp_identifier):
                dir_for_data[exp_identifier] = directory
            elif  dir_for_data[exp_identifier] != directory:
                print ">>> PROBLEM <<<<<same data set in different dirs ?"
                print exp_identifier, dir_for_data[exp_identifier], directory
    
    count = 0 # check on the tcga webpage how many
    for directory in newest_revision.keys():
        fields = directory.split('/')
        tumor_short = fields[0].upper()
        if tumor_short == 'READ': tumor_short = 'REA'

        path = local_dir + "/" + tumor_short
        check_and_make(path)
        path +=  "/CNV_SNP_Array"
        check_and_make(path)

        print directory, tumor_short, path
        for exp_identifier, revision in newest_revision[directory].iteritems():
            unzipped = ".".join([exp_identifier, str(revision), '0'])
            if  os.path.exists( path + "/" + unzipped):
                print "found", unzipped
                continue                
            filename = ".".join([exp_identifier, str(revision), '0', 'tar', 'gz' ])
            if  os.path.exists( path + "/" + filename):
                print "found", filename, path
                cmd = "tar -zxf " +  path +  "/" + filename + " -C " + path 
                print cmd
                retval = subprocess.call(cmd, shell=True)
                if retval==0: os.remove ( path +  "/" + filename)
                continue                
            
            full_url = root + "/" + directory + "/" + filename
            count += 1
            print count, path + "/" + filename
            # here comes the downloading part
            dwnldfile = urllib.URLopener()
            dwnldfile.retrieve(full_url, path + "/" + filename)

    print
    print "total: ", count


#########################################
if __name__ == '__main__':
    main()
