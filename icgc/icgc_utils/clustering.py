###############################
def parse_clust_out(clustering_output):
	isolated = []
	clusters = []

	inf = open(clustering_output,"r")
	reading_isolated = False
	reading_cluster = False
	list_of_res = []
	for line in inf:
		line = line.strip()
		if len(line)==0: continue
		if 'z-score' in line:
			zscore = float(line.split()[-1])
			if reading_cluster:
				clusters.append(list_of_res)
		elif 'isolated' in line:
			reading_isolated = True
			reading_cluster = False
			list_of_res = []
			continue
		elif  'cluster size' in line:
			if reading_isolated:
				isolated = list_of_res
			elif reading_cluster:
				clusters.append(list_of_res)
			reading_isolated = False
			reading_cluster = True
			list_of_res = []
			continue
		elif reading_cluster or reading_isolated:
			list_of_res.append(line)

	clusters = sorted(clusters, key=lambda cluster: -len(cluster))

	return zscore, isolated, clusters

