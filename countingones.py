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
PARA_BIAS       = 0.75
VIRULENCE       = 0.75
MUTATE_CHANCE   = 0.03
GENERATIONS     = 600
USE_SEED        = 3534566918815802929

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

def mainLoop(nIters, hostList, paraList, virulence):
  global relFitness
  for iteration in range(nIters):
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

# # max virulence run
# (hostList, paraList) = initLists()
# hostResultsList = []
# paraResultsList = []
# mainLoop(GENERATIONS, hostList, paraList, 1)
# plt.plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], 'o', color='red', markersize=0.1);
# plt.plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], 'o', color='blue', markersize=0.1);

#split virulence run
(hostList, paraList) = initLists()
hostResultsList = []
paraResultsList = []
relFitness = []
mainLoop(GENERATIONS, hostList, paraList, 1)
# mainLoop(GENERATIONS-200, hostList, paraList, VIRULENCE)
f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, gridspec_kw={'height_ratios': [5, 1, 1]})
ax1.xaxis.set_ticks_position('none')
ax1.plot([a[0] for a in hostResultsList], [a[1] for a in hostResultsList], 'o', color='red', markersize=0.1);
ax1.plot([a[0] for a in paraResultsList], [a[1] for a in paraResultsList], 'o', color='blue', markersize=0.1);
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