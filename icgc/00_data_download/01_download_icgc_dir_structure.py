#! /usr/bin/python

import urllib2

#########################################
def main():
	url = "https://dcc.icgc.org/api/v1/download/info?dir=releases&recursive=true"
	response = urllib2.urlopen(url)
	of = open ("dirstruct.json", "w")
	of.write(response.read())
	of.close()



#########################################
if __name__ == '__main__':
	main()
