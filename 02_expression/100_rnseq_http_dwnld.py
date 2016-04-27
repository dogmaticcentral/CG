#!/usr/bin/python
#
#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2016 Ivana Mihalek.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#
#
# run twice to unzip!
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

    file_names = open('expression_files.txt')
    all_files = {}
    for line in [x.rstrip().replace(' ','')  for x in file_names]:
        [directory, filename] = line.split(':')
        if directory[-1] == '/':   directory = directory[:-1]
        if not all_files.has_key(directory): all_files[directory] = []
        all_files[directory].append(filename)
    file_names.close()

    # TCGA has all revisions dumped in the same directory - see which one is the newest
    # here we have sever subdirs, correponding to different revisions
    revisions = {}
    for directory  in all_files.keys():
        fields = directory.split ('.')
        root_name = '.'.join(fields[:-2])
        revision  = fields[-2]
        if not revisions.has_key(root_name): revisions[root_name] = []
        revisions[root_name].append(int(revision))

    count = 0

    for root_name, revs in revisions.iteritems():
        latest_revision_dir = root_name + '.' + str(max(revs)) + '.0'

        fields = latest_revision_dir.split('/')
        tumor_short = fields[0].upper()
        if tumor_short == 'READ': tumor_short = 'REA'

        path = local_dir + "/" + tumor_short
        check_and_make(path)
        path += '/Expression_Genes'
        check_and_make(path)

        for file in all_files[latest_revision_dir]:

            full_url   = root + "/" + latest_revision_dir + "/" + file
            local_path = path + "/" + latest_revision_dir.split('/')[-1]
            check_and_make(local_path)

            # check if the file exists by any chance
            if  os.path.exists(local_path+ "/" + file):
                print file, "found in", local_path
            else:
                print "downloading " + file+ " to " + local_path
                dwnldfile = urllib.URLopener()
                dwnldfile.retrieve(full_url, local_path+ "/" + file)
            count += 1

    print
    print "total: ", count


#########################################
if __name__ == '__main__':
    main()
