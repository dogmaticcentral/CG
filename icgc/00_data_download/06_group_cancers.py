#! /usr/bin/python3
#
# This source code is part of icgc, an ICGC processing pipeline.
# 
# Icgc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Icgc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#

import os, subprocess

#########################################
from config import Config


def main():

	groups = {'PACA':['PAAD','PAEN','PACA'],
				'KIRC':['KIRC','CCSK'],
				'AML':['LAML','AML'],
				'COCA':['COAD','COCA','READ'],
				'ESAD':['ESAD','ESCA'],
				'LICA':['LIAD','LICA','LIRI','LINC','LIHC','LIHM'],
				'GACA':['STAD','GACA'],
				'MELA':['MELA','SKCM'],
				'PBCA':['PBCA', 'PEME'] }
	for group_head, group in groups.items():
		tmpdir = "{}.tmp".format(group_head)
		os.mkdir("{}/{}".format(Config.data_home_local, tmpdir))
		for directory in group:
			os.rename("{}/{}".format(Config.data_home_local, directory), "{}/{}/{}".format(Config.data_home_local, tmpdir, directory))
		os.rename("{}/{}".format(Config.data_home_local, tmpdir), "{}/{}".format(Config.data_home_local, group_head))

#########################################
if __name__ == '__main__':
	main()
