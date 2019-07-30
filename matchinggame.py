from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import copy
import sys
import pylab

BIT_LENGTH        = 100
POP_SIZE          = 50
HOST_BIAS         = 0.5
PARA_BIAS         = 0.5
MUTATE_CHANCE     = 0.03
GENERATIONS       = 600
USE_SEED          = None
SELECTION_SIZE    = 5
COMPETITION_SIZE  = 10
THRESHOLD         = 30
VIRULENCE         = 0.5
RESOURCE_SHARING  = True

class Participant():
  """Host or Parasite parent"""
  bitList = []
  score = 0.0
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
      host.victorIndices = []
      host.loserIndices = []
    for para in paraList:
      para.mutate()
      para.score = 0
      para.nLosses = 0

    #shuffle both lists and have the two populations compete pairwise
    hostWins = 0
    for x in range(COMPETITION_SIZE):
      shuffle(hostList)
      for popIndex in range(POP_SIZE):
        x  = paraList[popIndex].bitList.count(True)
        p1 = 0.5*(1.0+np.tanh((x-50)/7))
        matching = True if random()<p1 else False
        nMatching = 0
        for n in range(BIT_LENGTH):
          if paraList[popIndex].bitList[n] == matching and hostList[popIndex].bitList[n] == matching:
            nMatching += 1

        if nMatching >= THRESHOLD:
          #host wins
          if RESOURCE_SHARING:
            hostList[popIndex].loserIndices += [popIndex]
            paraList[popIndex].nLosses += 1
          else:
            hostList[popIndex].score += 1
          hostWins += 1
        else:
          #para wins
          if RESOURCE_SHARING:
            hostList[popIndex].victorIndices += [popIndex]
          else:
            paraList[popIndex].score += 1

    if RESOURCE_SHARING:
      for popIndex in range(POP_SIZE):
        #add score to paras
        for victorIndex in hostList[popIndex].victorIndices:
          paraList[victorIndex].score += 1.0/float(len(hostList[popIndex].victorIndices))
        #add score to hosts
        for loserIndex in hostList[popIndex].loserIndices:
          hostList[loserIndex].score += 1.0/paraList[loserIndex].nLosses

    relFitness.append([iteration,hostWins/(COMPETITION_SIZE*POP_SIZE)])

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
      hostSample = [hostList[i] for i in choices(range(POP_SIZE),k=SELECTION_SIZE)] #use choices for replacement
      best = hostSample[0]
      for host in hostSample:
        best = host if host.fitness > best.fitness else best
      newHostList.append(copy.deepcopy(best))

      paraSample = [paraList[i] for i in choices(range(POP_SIZE),k=SELECTION_SIZE)]
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
    hostList += [Participant(HOST_BIAS)]
    paraList += [Participant(PARA_BIAS)]
  return (hostList, paraList)


if USE_SEED is None:
  USE_SEED = randrange(sys.maxsize)
  rng = Random(USE_SEED)
print("Seed was:", USE_SEED)



fig = plt.figure(figsize=[10,10])
gs0 = fig.add_gridspec( 2, 2, width_ratios  = [1,1], height_ratios = [1,1], hspace=0.2, wspace=0.1)

gs = []
ax = []

for x in range(4):
  gs.append(gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[x], hspace=0.075, height_ratios=[5,1]))
  ax.append(plt.subplot(gs[x][0], xlim=[0,GENERATIONS], ylim=[0,BIT_LENGTH]))
  ax.append(plt.subplot(gs[x][1], xlim=[0,GENERATIONS], ylim=[0,1]))


#max virulence, no resource sharing
seed(USE_SEED)
VIRULENCE = 1.0
RESOURCE_SHARING = False
(hostList, paraList) = initLists()
hostResultsList = []
paraResultsList = []
relFitness = []
(hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, VIRULENCE)

ax[0].set_title("max virulence, no resource sharing".title(), fontsize=14)
ax[0].xaxis.set_visible(False)
ax[0].plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], 'o', color='red', markersize=0.1)
ax[0].plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], 'o', color='blue', markersize=0.1)
ax[0].axhline(y=PARA_BIAS*BIT_LENGTH, linestyle=':', color='black')

ax[1].tick_params(direction='in', left=True, right=True)
ax[1].plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)


#max virulence, with resource sharing
seed(USE_SEED)
VIRULENCE = 1.0
RESOURCE_SHARING = True
(hostList, paraList) = initLists()
hostResultsList = []
paraResultsList = []
relFitness = []
(hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, VIRULENCE)

ax[2].set_title("max virulence, with resource sharing".title(), fontsize=14)
ax[2].xaxis.set_visible(False)
ax[2].plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], 'o', color='red', markersize=0.1)
ax[2].plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], 'o', color='blue', markersize=0.1)
ax[2].axhline(y=PARA_BIAS*BIT_LENGTH, linestyle=':', color='black')

ax[3].tick_params(direction='in', left=True, right=True)
ax[3].plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)


#null virulence, no resource sharing
seed(USE_SEED)
VIRULENCE = 0.5
RESOURCE_SHARING = False
(hostList, paraList) = initLists()
hostResultsList = []
paraResultsList = []
relFitness = []
(hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, VIRULENCE)

ax[4].set_title("null virulence, no resource sharing".title(), fontsize=14)
ax[4].xaxis.set_visible(False)
ax[4].plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], 'o', color='red', markersize=0.1)
ax[4].plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], 'o', color='blue', markersize=0.1)
ax[4].axhline(y=PARA_BIAS*BIT_LENGTH, linestyle=':', color='black')

ax[5].tick_params(direction='in', left=True, right=True)
ax[5].plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)


#null virulence, with resource sharing
seed(USE_SEED)
VIRULENCE = 0.5
RESOURCE_SHARING = True
(hostList, paraList) = initLists()
hostResultsList = []
paraResultsList = []
relFitness = []
(hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, VIRULENCE)

ax[6].set_title("null virulence, with resource sharing".title(), fontsize=14)
ax[6].xaxis.set_visible(False)
ax[6].plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], 'o', color='red', markersize=0.1)
ax[6].plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], 'o', color='blue', markersize=0.1)
ax[6].axhline(y=PARA_BIAS*BIT_LENGTH, linestyle=':', color='black')

ax[7].tick_params(direction='in', left=True, right=True)
ax[7].plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)

plt.show()