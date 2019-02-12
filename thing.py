from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import pylab


class Participant():
    """Host or Parasite parent"""
    bitList = [False] * 100
    score = 0
    fitness = 0
    bias = 0.5

    def __init__(self, bias):
        self.bitList = [False] * 100
        self.bias = bias

    def mutate(self):
        newBitList = []
        for bit in self.bitList:
            if random() < 0.03:
                newBitList += [True] if random()<self.bias else [False]
            else:
                newBitList += [bit]
        self.bitList = newBitList

def compete(part1, part2):
    part1.score = 0
    part2.score = 0
    if sum(part1.bitList) > sum(part2.bitList):
        part1.score += 1
    elif sum(part1.bitList) == sum(part2.bitList):
        part1.score += 0.5
        part2.score += 0.5
    else: 
        part2.score += 1

hostList = []
paraList = []
hostResultsList = []
paraResultsList = []
def mainLoop(nIters=1):
    global hostList, paraList
    #initialise lists of participants
    for x in range(POP_SIZE):
        hostList     += [Participant(HOST_BIAS)]
        paraList += [Participant(PARA_BIAS)]

    for iteration in range(nIters):
        #mutate all participants and reset scores
        for host in hostList:
            host.mutate()
            host.score = 0
        for para in paraList:
            para.mutate()
            para.score = 0

        #compete with 5 random opponents
        for host in hostList:
            for paraIndex in sample(range(POP_SIZE),5):
                para = paraList[paraIndex]
                if host.bitList.count(True) > para.bitList.count(True):
                    host.score += 1
                    #print("host")
                elif host.bitList.count(True) == para.bitList.count(True):
                    host.score += 0.5
                    para.score += 0.5
                    #print("draw")
                else: 
                    para.score += 1
                    #print("para")
            host.fitness = host.score

        maxScore = max([para.score for para in paraList])

        #normalise scores and calculate fitness of parasites
        for para in paraList:
            para.score = para.score/maxScore
            para.fitness = ((2 * para.score) / (VIRULENCE)) + ((para.score * para.score) / (VIRULENCE * VIRULENCE))

        for host in hostList:
            hostResultsList.append([iteration,host.bitList.count(True)])
        for para in paraList:
            paraResultsList.append([iteration,para.bitList.count(True)])

        #asexual breeeding with tournement size 5
        newHostList = []
        newParaList = []
        for index in range(POP_SIZE):
            hostSample = [hostList[i] for i in sample(range(POP_SIZE),5)]
            best = hostSample[0]
            for host in hostSample:
                best = host if host.fitness > best.fitness else best
            newHostList.append(copy.deepcopy(best))

            paraSample = [paraList[i] for i in sample(range(POP_SIZE),5)]
            best = paraSample[0]
            for para in paraSample:
                best = para if para.fitness > best.fitness else best
            newParaList.append(copy.deepcopy(best))
        hostList = newHostList
        paraList = newParaList

def maxFitness(listParticipants):
    best = list(listParticipants)[0]
    for part in list(listParticipants):
        best = part if part.fitness > best.fitness else best

    return best


POP_SIZE            = 25
PARA_BIAS           = 0.8
HOST_BIAS           = 0.5
VIRULENCE           = .5
#randseed = randrange(sys.maxsize)
#rng = Random(randseed)
#print("Seed was:", randseed)
seed(8885448394944156554)

mainLoop(1000)
plt.plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], '.', color='red');
plt.plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], '.', color='blue');
plt.show()