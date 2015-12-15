#!/usr/bin/python

# 9th Module - Intro to Bioinformatics
# Bruno Azenha Goncalves
# ICMC - USP

NUMBER_OF_CLUSTERS = 2
ROWS = 2
COLS = 1

from numpy import array
from scipy.cluster.vq import vq, kmeans, whiten
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams


def getData(file_path):
	
	f = open(file_path, 'r')
	first_line = True	
	columns = 0
	leukemia_class = []
	features = []

	for line in f:
		entry = line.split(",")
		if first_line:
			columns = len(entry)
			first_line = False
			continue

		entry_list = []
		for i in range(columns):
			entry[i].strip()
			if i == columns-1:
				leukemia_class.append(int(entry[i].strip()))

			else:
				entry_list.append(float(entry[i]))

		features.append(entry_list)

	return leukemia_class, array(features)

def plotGeneExepression(features, idx, leukemia_class):
	matplotlib.rcParams['figure.figsize'] = 20, 15
	matplotlib.rcParams.update({'font.size': 8})

	for cluster in range(NUMBER_OF_CLUSTERS):
		

		false_count = 0
		plt.subplot(ROWS, COLS, cluster+1)
		for i in range(len(features)):
			if idx[i] == cluster:
				plt.plot(features[i])
				if leukemia_class[i] != cluster:
					false_count += 1

		if cluster == 1:
			plt.title("Has Leukemia - {} False Positives.".format(false_count))
		else:
			plt.title("Does NOT have Leukemia - {} False Negatives.".format(false_count))
		

	matplotlib.pyplot.gcf()
	plt.savefig("leukemia.png")


def main():

	print("... Clustering ...")
	leukemia, features = getData("leukemia.txt")
	print(leukemia)

	# Normalizes the data for a better clustering
	whitened = whiten(features)

	# Does the actual clustering
	centroids, _ = kmeans(whitened, NUMBER_OF_CLUSTERS)
	print("Done!")

	# Assing each gene to it's correct cluster
	idx, _ = vq(whitened, centroids)

	print("... Ploting graph ...")
	plotGeneExepression(features, idx, leukemia)
	print("Done!")




main()