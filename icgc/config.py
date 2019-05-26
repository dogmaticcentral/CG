
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
import os


class Config:

	ref_assembly = "hg19"

	icgc_release  = "27"
	data_home_local = "/storage/databases/icgc/v" + icgc_release

	mysql_conf_file = "/home/ivana/.tcga_conf"
	ucsc_mysql_conf_file  = "/home/ivana/.ucsc_mysql_conf"

	tcga_icgc_table_correspondence = {
		"ACC_somatic_mutations" : None,
		"ALL_somatic_mutations" : "ALL_simple_somatic",
		"BLCA_somatic_mutations": "BLCA_simple_somatic",
		"BRCA_somatic_mutations": "BRCA_simple_somatic",
		"CESC_somatic_mutations": "CESC_simple_somatic",
		"CHOL_somatic_mutations": None,
		"COAD_somatic_mutations": "COCA_simple_somatic",
		"DLBC_somatic_mutations": "DLBC_simple_somatic",
		"ESCA_somatic_mutations": "ESAD_simple_somatic",
		"GBM_somatic_mutations" : "GBM_simple_somatic",
		"HNSC_somatic_mutations": "HNSC_simple_somatic",
		"KICH_somatic_mutations": "KICH_simple_somatic",
		"KIRC_somatic_mutations": "KIRC_simple_somatic",
		"KIRP_somatic_mutations": "KIRP_simple_somatic",
		"LAML_somatic_mutations": "AML_simple_somatic",
		"LGG_somatic_mutations" : "LGG_simple_somatic",
		"LIHC_somatic_mutations": "LICA_simple_somatic",
		"LUAD_somatic_mutations": "LUAD_simple_somatic",
		"LUSC_somatic_mutations": "LUSC_simple_somatic",
		"MESO_somatic_mutations": None,
		"OV_somatic_mutations"  : "OV_simple_somatic",
		"PAAD_somatic_mutations": "PACA_simple_somatic",
		"PCPG_somatic_mutations": None,
		"PRAD_somatic_mutations": "PRAD_simple_somatic",
		"READ_somatic_mutations": "COCA_simple_somatic",
		"SARC_somatic_mutations": "SARC_simple_somatic",
		"SKCM_somatic_mutations": "MELA_simple_somatic",
		"STAD_somatic_mutations": "GACA_simple_somatic",
		"TGCT_somatic_mutations": None,
		"THCA_somatic_mutations": "THCA_simple_somatic",
		"THYM_somatic_mutations": None,
		"UCEC_somatic_mutations": "UCEC_simple_somatic",
		"UCS_somatic_mutations" : "UTCA_simple_somatic",
		"UVM_somatic_mutations" : None
	}

	def __init__(self):
		path = os.path.dirname(os.path.abspath(__file__))
		pathdirs = path.split("/")
		while pathdirs[-1]!= 'icgc': pathdirs.pop()
		self.rootpath = "/".join(pathdirs)


	# rbf is a small C program that runs the simulation
	# to evaluate Fisher-like probabilities for bins of uneven size (i.e probaility of being chosen)
	rbf_dir = "c-utils/random_binsize_fisher"
	rbf_path_relative = "%s/rbf" % rbf_dir
	def rbf_path(self):
		full_rbf_path = "{}/{}".format(self.rootpath, self.rbf_path_relative)
		if not os.path.exists(full_rbf_path) or not os.access(full_rbf_path, os.X_OK):
			print("rbf executable not found as", full_rbf_path)
			print("Perhaps it needs to be compiled by chdir-ing to $ICGC_HOME/%s and typing 'make'." % self.rbf_dir)
			exit()
		return full_rbf_path

	# pc outputs clusters on the structure given the input
	# selection of residues and the cutoff neighboring distance.
	# It also provides a z-score for the nonrandomness of the clustering
	# compared to random selection of the same size.
	pc_dir = "c-utils/clustering"
	pc_path_relative = "%s/pc" % pc_dir
	def pc_path(self):
		full_pc_path = "{}/{}".format(self.rootpath, self.pc_path_relative)
		if not os.path.exists(full_pc_path) or not os.access(full_pc_path, os.X_OK):
			print("pc executable not found as", full_pc_path)
			print("Perhaps it needs to be compiled by chdir-ing to $ICGC_HOME/%s and typing 'make'." % self.pc_dir)
			exit()
		return full_pc_path

