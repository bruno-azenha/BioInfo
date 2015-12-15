#!/usr/bin/python

# 9th Module - Intro to Bioinformatics
# Bruno Azenha Goncalves
# ICMC - USP

NUMBER_OF_CLUSTERS = 10
ROWS = 4
COLS = 3

X_AXIS = [40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260]

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
	genes = []
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
			if i == 0:
				genes.append(entry[i])

			else:
				entry_list.append(float(entry[i]))

		features.append(entry_list)

	return genes, array(features)

def plotGeneExepression(features, idx):
	matplotlib.rcParams['figure.figsize'] = 20, 15
	matplotlib.rcParams.update({'font.size': 8})

	for cluster in range(NUMBER_OF_CLUSTERS):
		count = 0
		plt.subplot(ROWS, COLS, cluster+1)
		plt.xlim((40,260))
		for i in range(len(features)):
			if idx[i] == cluster:
				plt.plot(X_AXIS, features[i])
				count += 1
		#plt.gca().axes.get_xaxis().set_visible(False)
		plt.title("{} Samples.".format(count))

	matplotlib.pyplot.gcf()
	plt.savefig("10Cluster.png")


def main():

	print("... Clustering ...")
	genes, features = getData("Spellman.csv")

	# Normalizes the data for a better clustering
	whitened = whiten(features)

	# Does the actual clustering
	centroids, _ = kmeans(whitened, NUMBER_OF_CLUSTERS)
	print("Done!")

	# Assing each gene to it's correct cluster
	idx, _ = vq(whitened, centroids)

	print("... Ploting graph ...")
	plotGeneExepression(features, idx)
	print("Done!")




main()