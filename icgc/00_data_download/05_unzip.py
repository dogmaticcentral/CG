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

	for root, dirs, files in os.walk(Config.data_home_local):
		for file in files:
			if file.endswith(".gz"):
				print(os.path.join(root, file))
				cmd = "gunzip " + os.path.join(root, file)
				subprocess.call(["bash","-c", cmd])

#########################################
if __name__ == '__main__':
	main()
