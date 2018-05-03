#! /usr/bin/python


import os, subprocess

#########################################
def main():
	data_home_local = "/data/icgc"
	for root, dirs, files in os.walk(data_home_local):
		for file in files:
			if file.endswith(".gz"):
				print os.path.join(root, file)
				cmd = "gunzip " + os.path.join(root, file)
				subprocess.call(["bash","-c", cmd])

#########################################
if __name__ == '__main__':
	main()
