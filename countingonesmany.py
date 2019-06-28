from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import pylab

BIT_LENGTH      = 100
POP_SIZE        = 25
HOST_BIAS       = 0.5
PARA_BIAS       = 0.9
VIRULENCE       = 1
MUTATE_CHANCE   = 0.03
GENERATIONS     = 600
USE_SEED        = None

class Participant():
  """Host or Parasite parent"""
  bitList = [False] * BIT_LENGTH
  score = 0
  fitness = 0
  bias = 0.5

  def __init__(self, bias):
    self.bitList = [False] * BIT_LENGTH
    self.bias = bias

  def mutate(self):
    newBitList = []
    for bit in self.bitList:
      if random() < MUTATE_CHANCE:
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

def mainLoop(nIters, iterOffset, hostList, paraList, virulence):
  global dGens
  for iteration in range(iterOffset, nIters):
    #mutate all participants and reset scores
    for host in hostList:
      host.mutate()
      host.score = 0
    for para in paraList:
      para.mutate()
      para.score = 0

    #shuffle both lists and have the two populations compete pairwise
    hostWins = 0
    for x in range(5):
      shuffle(hostList)
      for popIndex in range(POP_SIZE):
        if sum(hostList[popIndex].bitList) > sum(paraList[popIndex].bitList):
          hostList[popIndex].score += 1
          hostWins+=1
        elif sum(hostList[popIndex].bitList) == sum(paraList[popIndex].bitList):
          hostList[popIndex].score += 0.5
          paraList[popIndex].score += 0.5
          hostWins +=0.5
        else: 
          paraList[popIndex].score += 1
    dGens += 1 if (hostWins/(5*POP_SIZE)<0.025) else 0
    #normalise scores and calculate fitness of parasites
    maxScore = max([para.score for para in paraList])
    if maxScore > 0:
      for para in paraList:
        para.score = float(para.score)/maxScore
        para.fitness = ((2.0 * para.score) / (virulence)) - ((para.score * para.score) / (virulence * virulence))


    for host in hostList:
      host.fitness = host.score
      #print(host.score)
      hostResultsList.append([iteration,host.bitList.count(True)])
    for para in paraList:
      paraResultsList.append([iteration,para.bitList.count(True)])

    #asexual breeeding with tournement size 5
    newHostList = []
    newParaList = []
    for index in range(POP_SIZE):
      hostSample = [hostList[i] for i in choices(range(POP_SIZE),cum_weights=None,k=5)] #use choices for replacement
      best = hostSample[0]
      for host in hostSample:
        best = host if host.fitness > best.fitness else best
      newHostList.append(copy.deepcopy(best))

      paraSample = [paraList[i] for i in choices(range(POP_SIZE),cum_weights=None,k=5)]
      best = paraSample[0]
      for para in paraSample:
        best = para if para.fitness > best.fitness else best
      newParaList.append(copy.deepcopy(best))
    hostList = newHostList
    paraList = newParaList
  return (hostList, paraList)

def maxFitness(listParticipants):
  best = list(listParticipants)[0]
  for part in list(listParticipants):
    best = part if part.fitness > best.fitness else best

  return best
def initLists():  
  #initialise lists of participants
  hostList = []
  paraList = []
  for x in range(POP_SIZE):
    hostList   += [Participant(HOST_BIAS)]
    paraList += [Participant(PARA_BIAS)]
  return (hostList, paraList)


if USE_SEED is None:
  USE_SEED = randrange(sys.maxsize)
  rng = Random(USE_SEED)
print("Seed was:", USE_SEED)
seed(USE_SEED)

MUTATE_CHANCE = float(sys.argv[1])
VIRULENCE = float(sys.argv[2])

file = open("results_" + sys.argv[1] + "_" + sys.argv[2] + ".txt", "w+")

PARA_BIAS = 0.45
samplesize = 10
avgDGens = []
for y in range(11):
  dGens = 0
  PARA_BIAS += 0.04 if abs(PARA_BIAS-0.95)<0.01 else 0.05
  for x in range(samplesize):
    hostResultsList = []
    paraResultsList = []
    (hostList, paraList) = initLists()
    (hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, VIRULENCE)
  avgDGens.append([PARA_BIAS,dGens/samplesize]);
file.write(str(avgDGens))

plt.show()