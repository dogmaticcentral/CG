#! /usr/bin/python3

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
