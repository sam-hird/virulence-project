from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import copy
import sys
import pylab

BIT_LENGTH      = 100
POP_SIZE        = 25
HOST_BIAS       = 0.5
PARA_BIAS       = 0.75
VIRULENCE       = 1.0
MUTATE_CHANCE   = 0.03
GENERATIONS     = 600
USE_SEED        = 3534566918815802929
FIG             = 2

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


def mainLoop(nIters, iterOffset, hostList, paraList, virulence):
  global relFitness
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
    relFitness.append([iteration,hostWins/(5*POP_SIZE)])
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

if FIG == 1:
  PARA_BIAS = 0.75
elif FIG == 2 :
  PARA_BIAS = 0.9

#max virulence run
(hostList, paraList) = initLists()
hostResultsList = []
paraResultsList = []
relFitness = []
(hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, VIRULENCE)
fig = plt.figure(figsize=[15,5])
gs = fig.add_gridspec( 3, 3, width_ratios  = [1,1,1], height_ratios = [5,1,1])

ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[3], sharex=ax1)
ax3 = plt.subplot(gs[6], sharex=ax1)

ax1.set_title('Maximum Virulence', fontsize=16)
ax1.set_xlim([0,GENERATIONS])
ax1.set_ylim([0,BIT_LENGTH])
ax1.set_yticks(np.arange(0,BIT_LENGTH+1,step=25))
ax1.set_ylabel("Total number of 1s", fontsize=12)
ax1.xaxis.set_visible(False)
ax1.axhline(y=PARA_BIAS*BIT_LENGTH, linestyle=':', color='black')
ax1.axhline(y=HOST_BIAS*BIT_LENGTH, linestyle=':', color='black')
ax1.plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], 'o', color='red', markersize=0.1);
ax1.plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], 'o', color='blue', markersize=0.1);

ax2.set_ylim([0,1])
ax2.spines['bottom'].set_visible(False)
ax2.xaxis.set_visible(False)
ax2.tick_params(direction='in', left=True, right=True)
ax2.plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)

ax3.set_ylabel('               Relative fitness \n ', fontsize=12)
ax3.set_xlabel('Generations', fontsize=12)
ax3.set_ylim([0,1])
ax3.spines['top'].set_visible(False)
ax3.tick_params(direction='in', left=True, right=True)
ax3.plot([a[0] for a in relFitness], [1-a[1] for a in relFitness], '.', color="black", markersize=1)

#split virulence run
seed(USE_SEED)
(hostList, paraList) = initLists()
hostResultsList = []
paraResultsList = []
relFitness = []
(hostList, paraList) = mainLoop(350, 0, hostList, paraList, VIRULENCE)
VIRULENCE = 0.75
(hostList, paraList) = mainLoop(GENERATIONS, 350, hostList, paraList, VIRULENCE)

ax4 = plt.subplot(gs[1], sharey=ax1)
ax5 = plt.subplot(gs[4], sharex=ax4)
ax6 = plt.subplot(gs[7], sharex=ax4)

ax4.set_title('Maximum Virulence Gen 0-350 \n Moderate Virulence Gen 350-600', fontsize=16)
ax4.set_xlim([0,GENERATIONS])
ax4.set_yticks(np.arange(0,BIT_LENGTH+1,step=25))
ax4.xaxis.set_visible(False)
ax4.axhline(y=PARA_BIAS*BIT_LENGTH, linestyle=':', color='black')
ax4.axhline(y=HOST_BIAS*BIT_LENGTH, linestyle=':', color='black')
ax4.axvline(x=350, linestyle=':', color='black')
ax4.plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], 'o', color='red', markersize=0.1);
ax4.plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], 'o', color='blue', markersize=0.1);

ax5.set_ylim([0,1])
ax5.spines['bottom'].set_visible(False)
ax5.xaxis.set_visible(False)
ax5.tick_params(direction='in', left=True, right=True)
ax5.axvline(x=350, linestyle=':', color='black')
ax5.plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)

ax6.set_ylim([0,1])
ax6.spines['top'].set_visible(False)
ax6.tick_params(direction='in', left=True, right=True)
ax6.set_xlabel('Generations', fontsize=12)
ax6.axvline(x=350, linestyle=':', color='black')
ax6.plot([a[0] for a in relFitness], [1-a[1] for a in relFitness], '.', color="black", markersize=1)

#moderate virulence run
seed(USE_SEED)
(hostList, paraList) = initLists()
hostResultsList = []
paraResultsList = []
relFitness = []
(hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, VIRULENCE)

ax7 = plt.subplot(gs[2], sharey=ax1)
ax8 = plt.subplot(gs[5], sharex=ax7)
ax9 = plt.subplot(gs[8], sharex=ax7)

ax7.set_title('Moderate Virulence', fontsize=16)
ax7.set_xlim([0,GENERATIONS])
ax7.set_yticks(np.arange(0,BIT_LENGTH+1,step=25))
ax7.xaxis.set_visible(False)
ax7.axhline(y=PARA_BIAS*BIT_LENGTH, linestyle=':', color='black')
ax7.axhline(y=HOST_BIAS*BIT_LENGTH, linestyle=':', color='black')
ax7.plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], 'o', color='red', markersize=0.1);
ax7.plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], 'o', color='blue', markersize=0.1);

ax8.set_ylim([0,1])
ax8.spines['bottom'].set_visible(False)
ax8.xaxis.set_visible(False)
ax8.tick_params(direction='in', left=True, right=True)
ax8.plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)

ax9.set_ylim([0,1])
ax9.spines['top'].set_visible(False)
ax9.set_xlabel('Generations', fontsize=12)
ax9.tick_params(direction='in', left=True, right=True)
ax9.plot([a[0] for a in relFitness], [1-a[1] for a in relFitness], '.', color="black", markersize=1)

plt.show()