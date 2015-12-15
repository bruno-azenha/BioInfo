#!/usr/bin/python

# 6th Module - Intro to Bioinformatics
# Bruno Azenha Goncalves
# ICMC - USP
# Python program to build phylogenetic tree

from Bio import AlignIO
from Bio import Phylo
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

# Reads the alignment files
alignApe = AlignIO.read('genomes/human_primates_aligned.fasta', 'fasta')
alignHIV = AlignIO.read('genomes/alignedHIV.fasta', 'fasta')

# Creates the distance matrix
calculator = DistanceCalculator('ident')
dm_ape = calculator.get_distance(alignApe)
dm_hiv = calculator.get_distance(alignHIV)


# Jukes Cantor corrections
dm_ape_corrected = dm_ape
for d in dm_ape_corrected.matrix:
	d[:] = [-3/4*np.log(1-4/3*x) for x in d]

dm_hiv_corrected = dm_hiv
for d in dm_hiv_corrected.matrix:
	d[:] = [-3/4*np.log(1-4/3*x) for x in d]


# Constructs the tree using the upgma algorithm
constructor = DistanceTreeConstructor()

tree_ape = constructor.upgma(dm_ape)
tree_ape_corrected = constructor.upgma(dm_ape_corrected)

tree_hiv = constructor.upgma(dm_hiv)
tree_hiv_corrected = constructor.upgma(dm_hiv_corrected)

# Outputs the trees as a xml
Phylo.write(tree_ape, 'treeApe.xml', 'phyloxml')
Phylo.write(tree_ape_corrected, 'treeApe_corrected.xml', 'phyloxml')

Phylo.write(tree_hiv, 'treeHIV.xml', 'phyloxml')
Phylo.write(tree_hiv_corrected, 'treeHIV_corrected.xml', 'phyloxml')