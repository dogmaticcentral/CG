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
#A revision number can indicate the number of times an archive has been revised (starting from 0) and submitted to the DCC.
# However, the only requirement for revision numbers is that the revision number of the new archive is to be higher than
# that of the archive being replaced. Files that have been changed or added are captured in the changes and additions files, respectively
#
#  the series number should always be 0.
#
# note: run twice to get the unzipped files

import os, shutil


#########################################
def main():

    local_dir = '/Users/ivana/databases/TCGA'
    updating_somatic_mutations = True

    for tumor_short in os.listdir(local_dir):
        if 'usable' in tumor_short: continue
        print tumor_short
        path = local_dir + "/" + tumor_short + "/" + "Somatic_Mutations"
        if not os.path.exists(path): continue
        experiment_latest_revision = {}
        for data_set in os.listdir(path):
            # lets have this unzipped first
            if '.tar.gz' in data_set: continue
            fields = data_set.split ('.')
            revision = int(fields[-2])
            exp_identifier = '.'.join(fields[:-2])
            if not exp_identifier in experiment_latest_revision.keys():
                experiment_latest_revision[exp_identifier] = revision
            elif experiment_latest_revision[exp_identifier] < revision:
                experiment_latest_revision[exp_identifier] = revision
            # one more time:
        for data_set in os.listdir(path):
            if '.tar.gz' in data_set: continue

            fields = data_set.split ('.')
            revision = int(fields[-2])
            exp_identifier = '.'.join(fields[:-2])
            if revision == experiment_latest_revision[exp_identifier]: continue
            print data_set, "has newer revision:", experiment_latest_revision[exp_identifier]
            shutil.rmtree(path + "/" + data_set)

    # count = 0 # check on the tcga webpage how many
    # for directory in newest_revision.keys():
    #     print directory
    #     for exp_identifier, revision in newest_revision[directory].iteritems():
    #         print "\t", exp_identifier, revision
    #
    #     fields = directory.split('/')
    #     tumor_short = fields[0].upper()
    #     if tumor_short == 'READ': tumor_short = 'REA'
    #
    #     path = local_dir + "/" + tumor_short
    #     check_and_make(path)
    #     if updating_somatic_mutations:
    #         path += '/Somatic_Mutations'
    #     else:
    #         path +=  "/CNV_SNP_Array"
    #
    #     check_and_make(path)
    #
    #     for exp_identifier, revision in newest_revision[directory].iteritems():
    #         unzipped = ".".join([exp_identifier, str(revision), '0'])
    #         if  os.path.exists( path + "/" + unzipped):
    #             print "\tfound", unzipped
    #             continue
    #
    #         print "\tnew:  ", unzipped
    #         # unzip the new file
    #         filename = ".".join([exp_identifier, str(revision), '0', 'tar', 'gz' ])
    #         if  os.path.exists( path + "/" + filename):
    #             print "found", filename, path
    #             cmd = "tar -zxf " +  path +  "/" + filename + " -C " + path
    #             retval = subprocess.call(cmd, shell=True)
    #             if retval==0: os.remove ( path +  "/" + filename)
    #             continue
    #
    #         # if we got to here, the version of this file does not exist, zipped or unzipped
    #         # we are proceeding to download
    #         full_url = root + "/" + directory + "/" + filename
    #         count += 1
    #         print "downloading from ", full_url, " to ",  path + "/" + filename
    #         dwnldfile = urllib.URLopener()
    #         dwnldfile.retrieve(full_url, path + "/" + filename)
    #
    #         # is there an older revision of the same set of data by any chance?
    #         # lets move that to a new script, this is becoming to busy ...



#########################################
if __name__ == '__main__':
    main()
