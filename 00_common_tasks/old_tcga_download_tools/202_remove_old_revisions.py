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

    local_dir = '/mnt/databases/TCGA'
    updating_somatic_mutations = True

    for tumor_short in os.listdir(local_dir):
        if 'usable' in tumor_short: continue
        if 'auto' in tumor_short: continue
        print tumor_short
        path = local_dir + "/" + tumor_short + "/" + "Somatic_Mutations"
        if not os.path.exists(path): continue
        experiment_latest_revision = {}
        for data_set in os.listdir(path):
            # lets have this unzipped first
            if '.tar.gz' in data_set: continue
            fields = data_set.split ('.')
            print fields
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




#########################################
if __name__ == '__main__':
    main()
