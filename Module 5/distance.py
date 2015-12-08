#!/usr/bin/python

# 5th Module - Intro to Bioinformatics
# Bruno Azenha Goncalves
# ICMC - USP
# Little python program to calculate distances between aligned sequences

from Bio.Phylo.TreeConstruction import DistanceCalculator

# Reads a FASTA file, extracts and returns the information as string
def getGenome(filename):
    
    fastaFile = open(filename)

    # Throws the first line away
    fastaFile.readline()

    # Saves the rest to genome var stri
    genome = fastaFile.read()
    genome = genome.replace("\n", "")
    return genome

def main():
	
	calculator = DistanceCalculator()

	# Exercise 1
	print("Exercise 1:")
	genomeAfrican = getGenome("genomes/africanAligned.fasta")
	genomeIndian = getGenome("genomes/indianAligned.fasta")
	genomeMammoth = getGenome("genomes/mammothAligned.fasta")

	distAM = calculator._pairwise(genomeAfrican, genomeMammoth)
	distIM = calculator._pairwise(genomeIndian, genomeMammoth)

	print("Distance between African and Mammoth is {}.".format(distAM))
	print("Distance between Indian and Mammoth is {}.".format(distIM))

	# Exercise 3
	print("\nExercise 3:")
	genomeWhale = getGenome("genomes/whaleAligned.fasta")
	genomeCow = getGenome("genomes/cowAligned.fasta")
	genomeHippo = getGenome("genomes/hippoAligned.fasta")

	distWC = calculator._pairwise(genomeWhale, genomeCow)
	distWH = calculator._pairwise(genomeWhale, genomeHippo)
	
	print("Distance between Whale and Cow is {}.".format(distWC))
	print("Distance between Whale and Hippo is {}.".format(distWH))

main()