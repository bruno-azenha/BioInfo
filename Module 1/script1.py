#!/usr/bin/python

# 1st Module - Intro to Bioinformatics
# Bruno Azenha Goncalves
# ICMC - USP

import sys
import numpy as np
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

# Uses sliding window to reveal base composition
def slidingWindow(genome, windowSize):

    genomeLength = len(genome)
    if windowSize > genomeLength:
        print("Error: Window size is too big.")
        return 1
    
    # list that will store the GC frequency with the sliding window
    frequencyGC = []

    itNumber = 10
    for i in range(genomeLength-windowSize+1):
        counter = 0
        
        for j in range(0, windowSize):
            if genome[j+i] == 'C' or genome[j+i] == 'G':
                counter = counter + 1
        
        frequencyGC.append(counter)

    frequencyGC[:] = [float(x)/windowSize for x in frequencyGC]
    print(frequencyGC)
    return frequencyGC

def main():

    # Checks for the right command line argument number
    if len(sys.argv) != 3:
        print ("Usage: {0} [filepath] [windowSize]".format(sys.argv[0]))
        return 0

    genome = getGenome(sys.argv[1])
    windowSize = int(sys.argv[2])

    frequencyGC = slidingWindow(genome, windowSize)

    g = Gnuplot.Gnuplot(debug=1)
    g.title('Lambda GC content x Window Size({})'.format(windowSize))
    g('set yrange[0:1]')
    g('set term png')
    g('set output "labda{}.png"'.format(windowSize))
    g('set style data linespoints')
    g.plot(frequencyGC)
    
    return 0

main()