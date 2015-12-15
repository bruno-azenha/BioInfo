#!/usr/bin/python

# 6th Module - Intro to Bioinformatics
# Bruno Azenha Goncalves
# ICMC - USP

import sys
import numpy as np

# Reads a FASTA file, extracts and returns the information as string
def getGenome(filename):
    
    fastaFile = open(filename)

    # Throws the first line away
    fastaFile.readline()

    # Saves the rest to genome var stri
    genome = fastaFile.read()
    genome = genome.replace("\n", "")
    return genome

# Function that estimate Ka/Ks with pre aligned sequences
def NeiGojorobi(seq1, seq2):
	return (Ka/Ks)

def main():

	# ---- Exercise 1 ----
    print(" ------ Exercise 1 ------")
    print("")
    # String containing all the genome for each subject
    genomeHuman = getGenome("genomes/alignedHuman.fasta")
    genomeChimp = getGenome("genomes/alignedChimp.fasta")
    genomeMouse = getGenome("genomes/alignedMouse.fasta")

main()