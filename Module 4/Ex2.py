#!/usr/bin/python

# 4th Module - Intro to Bioinformatics
# Bruno Azenha Goncalves
# ICMC - USP

import sys
import numpy as np
import Gnuplot, Gnuplot.funcutils

WINDOW_SIZE = 20

# Returns a hidrophobiacity of a given aminoacid
def getHidrophobicity(aminoacid):
	a = aminoacid
	if a == "A":
		h = 1.8
	elif a == "R":
		h = -4.5
	elif a == "N":
		h = -3.5
	elif a == "D":
		h = -3.5
	elif a =="C":
		h = 2.5
	elif a == "Q":
		h = -3.5
	elif a == "E":
		h = -3.5
	elif a == "G":
		h = -0.4
	elif a == "H":
		h = -3.2
	elif a == "I":
		h = 1.5
	elif a == "L":
		h = 3.8
	elif a == "K":
		h = -3.9
	elif a == "M":
		h = 1.9
	elif a == "F":
		h = 2.8
	elif a == "P":
		h = -1.6
	elif a == "S":
		h = -0.8
	elif a == "T":
		h = -0.7
	elif a == "W":
		h = -0.9
	elif a == "Y":
		h = -1.3
	elif a == "V":
		h = 4.2
	else:
		exit("Error")
	return h

# Reads a FASTA file, extracts and returns the information as string
def getSequence(filename):
    
    fastaFile = open(filename)

    # Throws the first line away
    fastaFile.readline()

    # Saves the rest to genome var stri
    sequence = fastaFile.read()
    sequence = sequence.replace("\n", "")
    return sequence

def hydrophobicitySmooth(sequence):

	# First gets a list with each hidrophobicity
	hList = []
	for i in range(len(sequence)):
		hList.append(getHidrophobicity(sequence[i]))

	smoothList = []
	for i in range(len(hList) - WINDOW_SIZE):
		buff = 0;
		for j in range(WINDOW_SIZE):
			buff = buff + hList[i+j]

		smoothList.append(buff/WINDOW_SIZE)

	return smoothList


def main():

	sequence = getSequence("OR21.fasta")
	hsmooth = hydrophobicitySmooth(sequence)

	g = Gnuplot.Gnuplot(debug=1)
	g.title('Hidrophobicity x Position')
	g('set term png')
	g('set output "hydro.png"')
	g('set style data linespoints')
	g.plot(hsmooth)


main()
