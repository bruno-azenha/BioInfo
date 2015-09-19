#!/usr/bin/python

# 3rd Module - Intro to Bioinformatics
# Bruno Azenha Goncalves
# ICMC - USP

import sys
import numpy as np
from random import shuffle


# Reads a FASTA file, extracts and returns the information as string
def getGenome(filename):
    
    fastaFile = open(filename)

    # Throws the first line away
    fastaFile.readline()

    # Saves the rest to genome var stri
    genome = fastaFile.read()
    genome = genome.replace("\n", "")
    return genome

# Does the Needleman-Wunsch Algorithm
def Needle_Wunsch(seq1, seq2):

    MATCH = 5
    MISMATCH = -4
    GAP = -5

    # Initialize the table
    dpTable = np.zeros(shape = (len(seq1), len(seq2)), dtype = np.int16  )
    for i in range (len(seq1)):
        dpTable[i][0] = i * (-5)

    for j in range (len(seq2)):
        dpTable[0][j] = j * (-5)

    # Finds the Optimal Alignment
    for i in range (len(seq1))[1:]:
        for j in range (len(seq2))[1:]:    
                Match =  dpTable[i-1][j-1] + MATCH if seq1[i] == seq2[j] else dpTable[i-1][j-1] + MISMATCH
                Delete = dpTable[i-1][j] + GAP
                Insert = dpTable[i][j-1] + GAP
                dpTable[i][j] = max(Match, Insert, Delete)

    print(dpTable)
    
    # Now do the Traceback, starting from botton right
    alignedSeq1 = ""
    alignedSeq2 = ""

    i = len(seq1) - 1
    j = len(seq2) - 1

    while (i > 0 and j > 0):
        score = dpTable[i][j]
        scoreTop = dpTable[i][j-1]
        scoreLeft = dpTable[i-1][j]
        scoreDiag = dpTable[i-1][j-1]

        if score == scoreDiag + MATCH or score == scoreDiag + MISMATCH:
            alignedSeq1 = seq1[i] + alignedSeq1
            alignedSeq2 = seq2[j] + alignedSeq2
            i = i - 1
            j = j - 1

        elif score == scoreLeft + GAP:
            alignedSeq1 = seq1[i] + alignedSeq1
            alignedSeq2 = "-" + alignedSeq2
            i = i - 1

        elif score == scoreTop + GAP:
            alignedSeq1 = "-" + alignedSeq1
            alignedSeq2 = seq2[j] + alignedSeq2
            j = j - 1

        else:
            print("Wait, what?")

    # Complete the sequences:
    while (i >= 0):
        alignedSeq1 = seq1[i] + alignedSeq1
        alignedSeq2 = "-" + alignedSeq2
        i -= 1
    while (j >= 0):
        alignedSeq1 = "-" + alignedSeq1
        alignedSeq2 = seq2[j] + alignedSeq2
        j -= 1


    # Get statistics
    helper = ""
    matchNumber = 0
    gapNumber = 0
    mismatchNumber = 0

    for i in range(len(alignedSeq1)):
        if alignedSeq1[i] == alignedSeq2[i] and alignedSeq1[i] != "-":
            matchNumber = matchNumber + 1
            helper = helper + "|"
        elif alignedSeq1[i] == "-" or alignedSeq2[i] == "-":
            gapNumber = gapNumber + 1
            helper = helper + " "
        else:
            mismatchNumber = mismatchNumber + 1
            helper = helper + " "

    print("")
    print("Aniridia")
    print(alignedSeq1)
    print(helper)
    print(alignedSeq2)
    print("Drosophila")

    print("")
    print("Match Number: {}".format(matchNumber))
    print("Gap Number: {}".format(gapNumber))
    print("Mismatch Number: {}".format(mismatchNumber))

def main():

    genomeAniridia = getGenome("genomes/aniridia.fasta")
    genomeDrosophila = getGenome("genomes/drosophila.fasta")

    Needle_Wunsch(genomeAniridia, genomeDrosophila)


    return 1;

main()
