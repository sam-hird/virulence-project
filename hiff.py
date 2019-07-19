from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import pylab

BIT_LENGTH       = 64
POP_SIZE         = 25
HOST_BIAS        = 0.5
PARA_BIAS        = 0.5
MUTATE_CHANCE    = 0.03
GENERATIONS      = 600
USE_SEED         = None
VIRULENCE        = 0.5
SELECTION_SIZE   = 5
COMPETITION_SIZE = 5


class Participant():
  """Host or Parasite parent"""
  bitList = []
  score = 0.0
  fitness = 0
  bias = 0.5


  def __init__(self, bias):
    self.bitList = []
    for x in range(BIT_LENGTH):
      self.bitList += [1] if random()<bias else [0]
    self.bias = bias

  def mutate(self):
    newBitList = []
    for bit in self.bitList:
      if random() < MUTATE_CHANCE:
        newBitList += [1] if random()<self.bias else [0]
      else:
        newBitList += [bit]
    self.bitList = newBitList


def recursiveTransform(bitString):
  if len(bitString)==1:
    return bitString[0]
  else:
    sublistLength = len(bitString)//2
    return simpleTransform(recursiveTransform(bitString[0:sublistLength]), recursiveTransform(bitString[sublistLength:]))

def simpleTransform(leftBit, rightBit):
  if leftBit == None or rightBit == None or leftBit != rightBit:
    return None
  else:
    return leftBit

def recursiveFitness(bitString):
  if len(bitString) == 1:
    return simpleFitness(bitString[0])
  else:
    sublistLength = len(bitString)//2
    total =  len(bitString) * simpleFitness(recursiveTransform(bitString))
    total += recursiveFitness(bitString[:sublistLength])
    total += recursiveFitness(bitString[sublistLength:])
    return total

def simpleFitness(bit):
  return 1 if bit != None else 0

def maxfitness(maskedList):
  maxList = list(map(lambda x: 1 if x != None else None, maskedList))
  return recursiveFitness(maxList)

def maskLists(mask, bits):
  list1 = []
  list2 = []
  for x in range(len(mask)):
    if mask[x] == 1:
      list1 += [bits[x]]
      list2 += [None]
    else:
      list1 += [None]
      list2 += [bits[x]]
  return (list1, list2)


def mainLoop(nIters, iterOffset, hostList, paraList, virulence):
  global relFitness, connectednesses, maskLengths
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
    for x in range(COMPETITION_SIZE):
      shuffle(hostList)
      for popIndex in range(POP_SIZE):
        #mask lists
        (list1, list2) = maskLists(paraList[popIndex].bitList,hostList[popIndex].bitList)

        #record mask details
        connectedness = recursiveFitness(list(map(lambda x: 1 if x==1 else None, paraList[popIndex].bitList))) / recursiveFitness([1]*BIT_LENGTH)
        connectednesses.append([iteration,connectedness])
        maskLength = paraList[popIndex].bitList.count(1)
        maskLengths.append([iteration,maskLength])

        #get fitness proportions
        absHostHiffScore = recursiveFitness(list1)
        maxHostHiffScore = maxfitness(list1)
        hostScore = absHostHiffScore/maxHostHiffScore
        absHostHiffScores.append([iteration, absHostHiffScore])
        relHostHiffScores.append([iteration, hostScore])

        absParaHiffScore = recursiveFitness(list2)
        maxParaHiffScore = maxfitness(list2)
        paraScore = absParaHiffScore/maxParaHiffScore
        absParaHiffScores.append([iteration,absParaHiffScore])
        relParaHiffScores.append([iteration, paraScore])

        #see who wins
        if hostScore > paraScore:
          #host wins
          hostList[popIndex].score += 1
        elif paraScore > hostScore:
          #para wins

          paraList[popIndex].score += 1
        else:
          #draw
          hostList[popIndex].score += 0.5
          paraList[popIndex].score += 0.5
        

    relFitness.append([iteration,hostWins/(COMPETITION_SIZE*POP_SIZE)])

    #normalise scores and calculate fitness of parasites
    maxScore = max([para.score for para in paraList])
    if maxScore > 0:
      for para in paraList:
        para.score = float(para.score)/maxScore
        para.fitness = ((2.0 * para.score) / (virulence)) - ((para.score * para.score) / (virulence * virulence))


    for host in hostList:
      host.fitness = host.score

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
    hostList   += [Participant(HOST_BIAS)]
    paraList += [Participant(PARA_BIAS)]
  return (hostList, paraList)


if USE_SEED is None:
  USE_SEED = randrange(sys.maxsize)
  rng = Random(USE_SEED)
print("Seed was:", USE_SEED)
seed(USE_SEED)

hostList = []
paraList = []
maxScores           = [] #[[iteration, maxParaScore, maxHostScore]]   # maximum hiff score per generation per population
absHostHiffScores   = [] #[[iteration, HiffScore]]                    # actual hiff scores 
absParaHiffScores   = [] #[[iteration, HiffScore]]                    # actual hiff scores 
relHostHiffScores   = [] #[[iteration, HiffScore/maxScore]]           # percentage of possible hiff score for given mask
relParaHiffScores   = [] #[[iteration, HiffScore/maxScore]]           # percentage of possible hiff score for given mask
relFitness          = [] #[[iteration, relFitness]]                   # percentage of competitions won by host
connectednesses     = [] #[[iteration, connectedness]]
maskLengths         = [] #[[iteration, maskLength]]


(hostList, paraList) = initLists()
(hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, VIRULENCE)


f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, gridspec_kw={'height_ratios': [5, 5, 1]})
ax1.xaxis.set_ticks_position('none')
ax1.plot([a[0] for a in absHostHiffScores], [a[1] for a in absHostHiffScores], 'o', color='red', markersize=0.1);
ax1.plot([a[0] for a in absParaHiffScores], [a[1] for a in absParaHiffScores], 'o', color='blue', markersize=0.1);
ax2.plot([a[0] for a in connectednesses], [a[1] for a in connectednesses], 'o', color='red', markersize=0.1);
ax2.plot([a[0] for a in maskLengths], [a[1] for a in maskLengths], 'o', color='blue', markersize=0.1);

ax3.spines['top'].set_visible(False)
ax3.tick_params(direction='in', left=True, right=True)
ax3.plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)
plt.show()