from random import *
from operator import attrgetter

pop_size            = 25
para_bias           = 0.75
host_bias           = 0.5
virulence           = 0.75

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
    if sum(part1.bitList) > sum(part2.bitList):
        part1.score += 1
    elif sum(part1.bitList) == sum(part2.bitList):
        part1.score += 0.5
        part2.score += 0.5
    else: 
        part2.score += 1

hostList     = []
parasiteList = []
def mainLoop(nIters=1):
    global hostList, parasiteList
    #initialise lists of participants
    for x in range(pop_size):
        hostList     += [Participant(host_bias)]
        parasiteList += [Participant(para_bias)]

    for iteration in range(nIters):
        #mutate all participants
        for host in hostList + parasiteList:
            host.mutate()

        #compete with 5 random opponents
        for host in hostList:
            for paraIndex in sample(range(pop_size),5):
                compete(host,parasiteList[paraIndex])
            host.fitness = host.score                

        #calculate fitness of parasites
        for para in parasiteList:
            para.fitness = ((2 * para.score) / (virulence)) + ((para.score * para.score) / (virulence * virulence))

        #asexual breeeding with tournement size 5
        newHostList = []
        newParaList = []
        for index in range(pop_size):
            newHostList.append(maxFitness(hostList[i] for i in sample(range(pop_size),5)))
            newParaList.append(maxFitness(parasiteList[i] for i in sample(range(pop_size),5)))
        hostList     = newHostList
        parasiteList = newParaList

def maxFitness(listParticipants):
    best = Participant(0)
    for part in listParticipants:
        best = part if part.fitness > best.fitness else best
    return best