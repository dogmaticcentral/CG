#! /usr/bin/python3


import os, subprocess

#########################################
def main():
	release  = "27"
	data_home_local = "/storage/databases/icgc/v"+release
	for root, dirs, files in os.walk(data_home_local):
		for file in files:
			if file.endswith(".gz"):
				print(os.path.join(root, file))
				cmd = "gunzip " + os.path.join(root, file)
				subprocess.call(["bash","-c", cmd])

#########################################
if __name__ == '__main__':
	main()
