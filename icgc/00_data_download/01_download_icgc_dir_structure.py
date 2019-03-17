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
DEVNULL = open(os.devnull, 'wb')

# not sure how to download the listing of controlled files
# I am going with the assumption that for every 'open' file
# there is a 'controlled' version
#########################################
def main():

	url = "https://dcc.icgc.org/api/v1/download/info?dir=releases&recursive=true"

	#-L, --locationf the server reports that the requested page has moved to a different location
	# (indicated with a Location: header and a 3XX response code),
	# this option will make curl redo the request on the new place.
	cmd = "curl -L '{}' -o dirstruct.json ".format(url)
	subprocess.Popen(["bash", "-c", cmd])



#########################################
if __name__ == '__main__':
	main()
