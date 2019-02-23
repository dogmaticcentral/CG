#! /usr/bin/python3


import os, subprocess

#########################################
def main():
	release  = "27"
	data_home_local = "/storage/databases/icgc/v"+release
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
		os.mkdir("{}/{}".format(data_home_local, tmpdir))
		for directory in group:
			os.rename("{}/{}".format(data_home_local, directory), "{}/{}/{}".format(data_home_local, tmpdir, directory))
		os.rename("{}/{}".format(data_home_local, tmpdir), "{}/{}".format(data_home_local, group_head))

#########################################
if __name__ == '__main__':
	main()
