#! /usr/bin/python3


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
