#! /usr/bin/python

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
