#
# H-IFF style 2
# using mask and complement of mask for scoring
#


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
VIRULENCE        = 1
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

def maxFitness(maskedList):
  maxList = list(map(lambda x: 1 if x != None else None, maskedList))
  return recursiveFitness(maxList)

def maskLists(mask, bits):
  list1 = []
  list2 = []
  for x in range(len(mask)):
    if mask[x] == 0:
      list1 += [bits[x]]
      list2 += [None]
    else:
      list1 += [None]
      list2 += [bits[x]]
  return (list1, list2)


def mainLoop(nIters, iterOffset, hostList, paraList, virulence):
  global hostAbsHiffScores, hostRelHiffScores, hostNum1s
  global paraAbsHiffScores, paraRelHiffScores, paraNum1s

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
        connectedness = recursiveFitness(list(map(lambda x: 1 if x==1 else None, paraList[popIndex].bitList))) / recursiveFitness([1]*paraList[popIndex].bitList.count(1))
        connectednesses.append([iteration,connectedness])
        maskLength = paraList[popIndex].bitList.count(1)
        maskLengths.append([iteration,maskLength])

        #calculate scores
        hostAbsHiffScore = recursiveFitness(list1)
        hostScore = hostAbsHiffScore/maxFitness(list1)
        paraAbsHiffScore = recursiveFitness(list2)
        paraScore = paraAbsHiffScore/maxFitness(list2)

        #store results
        hostAbsHiffScores.append([iteration, hostAbsHiffScore])
        hostRelHiffScores.append([iteration, hostScore])
        paraAbsHiffScores.append([iteration, paraAbsHiffScore])
        paraRelHiffScores.append([iteration, paraScore])

        #see who wins
        if hostScore > paraScore:
          #host wins
          hostList[popIndex].score += 1
          hostWins +=1
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
      hostNum1s.append([iteration,host.bitList.count(1)])
    for para in paraList:
      paraNum1s.append([iteration,para.bitList.count(1)])

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

# hostList = []
# paraList = []
# maxScores           = [] #[[iteration, maxParaScore, maxHostScore]]   # maximum hiff score per generation per population

# hostNum1s           = [] #[[iteration, num1s]]
# paraNum1s           = [] #[[iteration, num1s]]
# hostAbsHiffScores   = [] #[[iteration, HiffScore]]                    # actual hiff scores 
# paraAbsHiffScores   = [] #[[iteration, HiffScore]]                    # actual hiff scores 
# hostRelHiffScores   = [] #[[iteration, HiffScore/maxScore]]           # percentage of possible hiff score for given mask
# paraRelHiffScores   = [] #[[iteration, HiffScore/maxScore]]           # percentage of possible hiff score for given mask
# relFitness          = [] #[[iteration, relFitness]]                   # percentage of competitions won by host
# connectednesses     = [] #[[iteration, connectedness]]
# maskLengths         = [] #[[iteration, maskLength]]

# #initialize lists of participants
# (hostList, paraList) = initLists()
# #do 600 iterations
# (hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, VIRULENCE)

# #define sublot parameters
# f, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(5, 1, sharex=True, subplot_kw={'xlim': {0,GENERATIONS}}, gridspec_kw={'height_ratios': [5, 5, 5, 5, 1]}, figsize=(5, 15))

# #plot num1s
# ax0.title.set_text("Proportion of 1s in genotype")
# ax0.set_ylim([0,1])
# ax0.plot([a[0] for a in hostNum1s], [a[1]/BIT_LENGTH for a in hostNum1s], 'o', color='red', markersize=0.1);
# ax0.plot([a[0] for a in paraNum1s], [a[1]/BIT_LENGTH for a in paraNum1s], 'o', color='blue', markersize=0.1);

# #plot abs hiff
# ax1.title.set_text("Absolute HIFF Scores")
# ax1.set_ylim([0,maxFitness([1]*BIT_LENGTH)])
# ax1.plot([a[0] for a in hostAbsHiffScores], [a[1] for a in hostAbsHiffScores], 'o', color='red', markersize=0.1);
# ax1.plot([a[0] for a in paraAbsHiffScores], [a[1] for a in paraAbsHiffScores], 'o', color='blue', markersize=0.1);

# #plot relative hiff
# ax2.title.set_text("Proportion of possible HIFF score")
# ax2.set_ylim([0,1])
# ax2.plot([a[0] for a in hostRelHiffScores], [a[1] for a in hostRelHiffScores], 'o', color='red', markersize=0.1, clip_on=False, zorder=100);
# ax2.plot([a[0] for a in paraRelHiffScores], [a[1] for a in paraRelHiffScores], 'o', color='blue', markersize=0.1, clip_on=False, zorder=100);

# #plot conectedness
# ax3.title.set_text("Connectedness of mask")
# ax3.set_ylim([0,1])
# ax3.plot([a[0] for a in connectednesses], [a[1] for a in connectednesses], 'o', color='blue', markersize=0.1);

# #plot relative fitness
# ax4.title.set_text("proportion of host victories")
# ax4.set_ylim([0,1])
# ax4.plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)

# #display graph
# plt.show()
