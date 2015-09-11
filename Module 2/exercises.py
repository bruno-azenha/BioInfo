#!/usr/bin/python

# 2nd Module - Intro to Bioinformatics
# Bruno Azenha Goncalves
# ICMC - USP

import sys
import numpy as np
from random import shuffle
import Gnuplot, Gnuplot.funcutils


# Reads a FASTA file, extracts and returns the information as string
def getGenome(filename):
    
    fastaFile = open(filename)

    # Throws the first line away
    fastaFile.readline()

    # Saves the rest to genome var stri
    genome = fastaFile.read()
    genome = genome.replace("\n", "")
    return genome

# Reverses the genome string
# - NOT USED -
def getReverse(genome):
    reverse = genome

    # "Buffers"
    reverse.replace("A", "a");
    reverse.replace("T", "t");
    reverse.replace("C", "c");
    reverse.replace("G", "g");

    # Reverses
    reverse.replace("a", "T")
    reverse.replace("t", "A")
    reverse.replace("c", "G")
    reverse.replace("g", "C")

    return reverse

# Searches for ORF
def SearchORF(genome, index):
    begin = index
    end = len(genome)

    orfList = [];
    numberOfBasesInORF = 0;
    maxOrfInCodons = 0

    # -1 for not found, index otherwise
    foundStart = -1

    for i in range(begin, end, 3):

        # Searches for start codon
        if (foundStart == -1 and genome[i:i+3] == "ATA"):
            #print("found start codon at {}".format(i))
            foundStart = i;

        # Searches for stop codon 
        if (foundStart != -1 and genome[i:i+3] in {"AGA", "AGG"}):
            #print("found stop codon at {}".format(i))
            orfList.append([foundStart, i])
            numberOfBasesInORF = numberOfBasesInORF + (i - foundStart) + 3
            if (maxOrfInCodons < (i - foundStart)/3 + 1):
                maxOrfInCodons = (i - foundStart)/3 + 1
            foundStart = -1

    return orfList, numberOfBasesInORF, maxOrfInCodons

def randomizeGenome(baseGenome):
    l = list(baseGenome)
    shuffle(l)
    randomizedGenome = ''.join(l)
    return randomizedGenome

def main():

    # ---- Exercise 1 ----
    print(" ------ Exercise 1 ------")
    print("")
    # String containing all the genome for each subject
    genomeHuman = getGenome("genomes/humanMT.fasta")
    genomeChimp = getGenome("genomes/chimpMT.fasta")
    genomeMouse = getGenome("genomes/mouseMT.fasta")


    # Human Mitochondria
    orfList, numberOfBasesInORF, maxORF = SearchORF(genomeHuman, 0);
    ratio = 100 * numberOfBasesInORF / float(len(genomeHuman))
    print("Human Mithocondria Genome:")
    print("{} ORFs were found.".format(len(orfList)))
    print("{}% is candidate protein.".format(ratio))
    print("")

    # Chimpanzee Mitochondria
    orfList, numberOfBasesInORF, maxORF = SearchORF(genomeChimp, 0);
    ratio = 100 * numberOfBasesInORF / float(len(genomeChimp))
    print("Chimpanzee Mithocondria Genome:")
    print("{} ORFs were found.".format(len(orfList)))
    print("{}% is candidate protein.".format(ratio))
    print("")

    # Mouse Mitochondria
    orfList, numberOfBasesInORF, maxORF = SearchORF(genomeMouse, 0);
    ratio = 100 * numberOfBasesInORF / float(len(genomeMouse))
    print("Mouse Mithocondria Genome:")
    print("{} ORFs were found.".format(len(orfList)))
    print("{}% is candidate protein.".format(ratio))
    print("")

    # ---- Exercise 2 ----
    print(" ------ Exercise 2 ------")
    print("")

    # Random Genome from Human Mitochondria
    genomeRandomHuman = randomizeGenome(genomeHuman)
    orfList, numberOfBasesInORF, maxORF = SearchORF(genomeRandomHuman, 0);
    ratio = 100 * numberOfBasesInORF / float(len(genomeRandomHuman))
    print("Randomly Generated Genome (Shuffling Human Mitochondria):")
    print("{} ORFs were found.".format(len(orfList)))
    print("{}\% is candidate protein.".format(ratio))
    print("Largest ORF is {} codons long.".format(maxORF))
    print("")

    # Random Genome from Chimpanzee Mitochondria
    genomeRandomChimp = randomizeGenome(genomeChimp)
    orfList, numberOfBasesInORF, maxORF = SearchORF(genomeRandomChimp, 0);
    ratio = 100 * numberOfBasesInORF / float(len(genomeRandomChimp))
    print("Randomly Generated Genome (Shuffling Chimpanzee Mitochondria):")
    print("{} ORFs were found.".format(len(orfList)))
    print("{}\% is candidate protein.".format(ratio))
    print("Largest ORF is {} codons long.".format(maxORF))
    print("")
    print("")

    # Random Genome from Mouse Mitochondria
    genomeRandomMouse = randomizeGenome(genomeMouse)
    orfList, numberOfBasesInORF, maxORF = SearchORF(genomeRandomMouse, 0);
    ratio = 100 * numberOfBasesInORF / float(len(genomeRandomMouse))
    print("Randomly Generated Genome (Shuffling Mouse Mitochondria):")
    print("{} ORFs were found.".format(len(orfList)))
    print("{}\% is candidate protein.".format(ratio))
    print("Largest ORF is {} codons long.".format(maxORF))
    print("")
    

    # ---- Exercise 3 ----

    genomeFlu = getGenome("genomes/fluGenome.fasta")
    genomeRandomFlu = randomizeGenome(genomeFlu)

    # Only considers ORFs longer than the LONGEST randomly generated ORF
    randomOrfList, numberOfBasesInORF, maxORF = SearchORF(genomeRandomFlu, 0)
    randomOrfList.sort(key=lambda randomOrfList: (randomOrfList[1]-randomOrfList[0])/3 + 2, reverse=True)  # sorts in place
    maxRandomOrf = (randomOrfList[0][1] - randomOrfList[0][0]) / 3 + 2

    orfList, numberOfBasesInORF, maxORF = SearchORF(genomeFlu, 0)
    ORFBefore = len(orfList)
    genomeBuffer = [];
    for i in range(len(orfList)):
        if ((orfList[i][1]- orfList[i][0])/3 + 2 > maxRandomOrf):
            genomeBuffer.append(orfList[i])

    print("")
    print ("Considers only ORFs longer than Longest random ORF")
    print ("Cutoff Codon Number: {}".format(maxRandomOrf))
    print ("ORFs before cutoff: {}".format(ORFBefore))
    print ("ORFs after cutoff: {}".format(len(genomeBuffer)))

    # Now consider ORFs longer than the 5% Longest randomly generated ORFs
    length = len(randomOrfList)
    pos = length/20 # Finds the top 5%
    maxRandomOrf = (randomOrfList[pos][1] - randomOrfList[pos][0]) / 3 + 2

    genomeBuffer = [];
    for i in range(len(orfList)):
        if ((orfList[i][1]- orfList[i][0])/3 + 2 > maxRandomOrf):
            genomeBuffer.append(orfList[i])

    print("")
    print ("Considers only ORFs longer than the 5% Longest random ORF")
    print ("Cutoff Codon Number: {}".format(maxRandomOrf))
    print ("ORFs before cutoff: {}".format(ORFBefore))
    print ("ORFs after cutoff: {}".format(len(genomeBuffer)))

     # Now consider ORFs longer than the 10% Longest randomly generated ORFs
    length = len(randomOrfList)
    pos = length/10 # Finds the top 10%
    maxRandomOrf = (randomOrfList[pos][1] - randomOrfList[pos][0]) / 3 + 2

    genomeBuffer = [];
    for i in range(len(orfList)):
        if ((orfList[i][1]- orfList[i][0])/3 + 2 > maxRandomOrf):
            genomeBuffer.append(orfList[i])

    print("")
    print ("Considers only ORFs longer than the 10% Longest random ORF")
    print ("Cutoff Codon Number: {}".format(maxRandomOrf))
    print ("ORFs before cutoff: {}".format(ORFBefore))
    print ("ORFs after cutoff: {}".format(len(genomeBuffer)))

main()


