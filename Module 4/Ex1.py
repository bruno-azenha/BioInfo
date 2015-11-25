#!/usr/bin/python

# 4th Module - Intro to Bioinformatics
# Bruno Azenha Goncalves
# ICMC - USP

RESULT_SEQUENCE = "01011101001"
HIDDEN_STATES =   "FFFBBBBBFFF"
NUMBER_OF_TRIALS = 10000

import sys
from random import randint

# Gets the result of a coin toss
# Can be a fair coin or not
def CoinToss(isFair):
	rand = randint(1,4)

	if isFair:
		coin = "0" if rand <= 2 else "1"
	else:
		coin = "0" if rand <= 3 else "1"
	
	return coin

# Returns True if the result is equal to the expected result sequence
# Returns False otherwise
def RunExperiment():

	# Initialize lists to store results
	coinList = ""

	# Throws a coin numberOfThrows times
	for i in range(len(HIDDEN_STATES)):
		coin = CoinToss(True if HIDDEN_STATES[i] == "F" else False)
		coinList = coinList + coin
		if coinList[i] != RESULT_SEQUENCE[i]: # Early discovery of mismatch
			return False

	return True

# Counts and prints the chance of the result sequence to occur given the hidden states.
def main():

	successes = 0

	for i in range (NUMBER_OF_TRIALS):
		if RunExperiment() == True:
			successes = successes + 1

	print ("Looks like the probability is {}.".format(float(successes)/NUMBER_OF_TRIALS))
	return 0

main()
