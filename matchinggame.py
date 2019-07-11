from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import pylab

BIT_LENGTH      = 100
POP_SIZE        = 50
HOST_BIAS       = 0.5
PARA_BIAS       = 0.5
MUTATE_CHANCE   = 0.003
GENERATIONS     = 600
USE_SEED        = None

class Participant():
  """Host or Parasite parent"""
  bitList = [False] * BIT_LENGTH
  score = 0
  fitness = 0
  bias = 0.5


  def __init__(self, bias):
    self.bitList = []
    for x in range(BIT_LENGTH):
      self.bitList += [False] if random()<bias else [True]
    self.bias = bias

  def mutate(self):
    newBitList = []
    for bit in self.bitList:
      if random() < MUTATE_CHANCE:
        newBitList += [True] if random()<self.bias else [False]
      else:
        newBitList += [bit]
    self.bitList = newBitList

def mainLoop(nIters, iterOffset, hostList, paraList, virulence):
  global relFitness
  for iteration in range(iterOffset, nIters):
    #mutate all participants and reset scores
    for host in hostList:
      host.mutate()
      host.score = 0
      host.victors = []
    for para in paraList:
      para.mutate()
      para.score = 0
      para.victors = []

    #shuffle both lists and have the two populations compete pairwise
    hostWins = 0
    for x in range(10):
      shuffle(hostList)
      for popIndex in range(POP_SIZE):
        x  = paraList[popIndex].bitList.count(True)
        p1 = 0.5*(1+np.tanh((x-50)/7))
        matching = True if random()<p1 else False
        nMatching = 0
        for n in range(BIT_LENGTH):
          if paraList[popIndex].bitList[n] == matching and hostList[popIndex].bitList[n] == matching:
            nMatching += 1

        if nMatching >= 30:
          #host wins
          if RESOURCE_SHARING:
            paraList[popIndex].victors += [hostList[popIndex]]
          else:
            hostList[popIndex].score += 1
          hostWins += 1
        else:
          #para wins
          if RESOURCE_SHARING:
            hostList[popIndex].victors += [paraList[popIndex]]
          else:
            paraList[popIndex].score += 1

    if RESOURCE_SHARING:
      for popIndex in range(POP_SIZE):
        #add score to paras
        for victor in hostList[popIndex].victors:
          victor.score += 1/len(hostList[popIndex].victors)
        #add score to hosts
        for victor in paraList[popIndex].victors:
          victor.score += 1/len(paraList[popIndex].victors)

    relFitness.append([iteration,hostWins/(10*POP_SIZE)])

    #normalise scores and calculate fitness of parasites
    maxScore = max([para.score for para in paraList])
    if maxScore > 0:
      for para in paraList:
        para.score = float(para.score)/maxScore
        para.fitness = ((2.0 * para.score) / (virulence)) - ((para.score * para.score) / (virulence * virulence))

    for host in hostList:
      host.fitness = host.score
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

def initLists():  
  #initialise lists of participants
  hostList = []
  paraList = []
  for x in range(POP_SIZE):
    hostList   += [Participant(HOST_BIAS)]
    paraList += [Participant(PARA_BIAS)]
  return (hostList, paraList)

USE_SEED = 4647813385336841684

if USE_SEED is None:
  USE_SEED = randrange(sys.maxsize)
  rng = Random(USE_SEED)
print("Seed was:", USE_SEED)
seed(USE_SEED)

VIRULENCE = 1
RESOURCE_SHARING = True

#split virulence run
(hostList, paraList) = initLists()
hostResultsList = []
paraResultsList = []
relFitness = []
(hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, VIRULENCE)
# (hostList, paraList) = mainLoop(GENERATIONS, 250, hostList, paraList, VIRULENCE)
f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, gridspec_kw={'height_ratios': [5, 1, 1]})
ax1.xaxis.set_ticks_position('none')
ax1.plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], 'o', color='red', markersize=0.1);
ax1.plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], 'o', color='blue', markersize=0.1);
ax1.set_ylim([0,100])
ax1.set_xlim([0,GENERATIONS])
ax2.set_ylim([0,1])
ax2.spines['bottom'].set_visible(False)
ax2.xaxis.set_ticks_position('none')
ax2.tick_params(direction='in', left=True, right=True)
ax2.plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)
ax3.set_ylim([0,1])
ax3.spines['top'].set_visible(False)
ax3.tick_params(direction='in', left=True, right=True)
ax3.plot([a[0] for a in relFitness], [1-a[1] for a in relFitness], '.', color="black", markersize=1)
plt.show()