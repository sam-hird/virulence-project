#
# H-IFF Style 1
# using a threshold for scoring
#

from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import copy
import sys
import pylab

BIT_LENGTH            = 32
POP_SIZE              = 25
HOST_BIAS             = 0.5
PARA_BIAS             = 0.05
HOST_MUTATE_CHANCE    = 0.02
PARA_MUTATE_CHANCE    = 0.02
GENERATIONS           = 1000
USE_SEED              = None
VIRULENCE             = 1
SELECTION_SIZE        = 10
COMPETITION_SIZE      = 10
THRESHOLD             = 1.0
INIT_METHOD           = 2


class Participant():
  #Host or parasite parent

  #default values for parameters
  bitList = []
  score = 0.0
  fitness = 0
  bias = 0.5
  mutationRate = 0.02


  def __init__(self, bias, mutationRate, bitList):
    if bitList == []:
      self.bitList = []
      for x in range(BIT_LENGTH):
        self.bitList += [1] if random()<bias else [0]
    else:
      self.bitList = bitList

    self.mutationRate = mutationRate
    self.bias = bias

  def mutate(self):
    newBitList = []
    for bit in self.bitList:
      if random() < self.mutationRate:
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
    if mask[x] == 1:
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
        maskHIFF = recursiveFitness(list(map(lambda x: 1 if x==1 else None, paraList[popIndex].bitList)))
        maskHIFFMax = 0 if paraList[popIndex].bitList.count(1) == 0 else recursiveFitness([1]*paraList[popIndex].bitList.count(1))
        connectedness =   0 if maskHIFFMax == 0 else maskHIFF / maskHIFFMax
        connectednesses.append([iteration,connectedness])

        #calculate scores
        hostAbsHiffScore = recursiveFitness(list1)
        hostScore =  0 if maxFitness(list1) == 0 else hostAbsHiffScore / maxFitness(list1)

        #store results
        hostAbsHiffScores.append([iteration, hostAbsHiffScore])
        hostRelHiffScores.append([iteration, hostScore])

        #see who wins
        if hostScore >= THRESHOLD:
          #host wins
          hostList[popIndex].score += 1
          hostWins +=1
        else:
          #para wins
          paraList[popIndex].score += 1

    for host in hostList:
      host.fitness = host.score
      hostNum1s.append([iteration,host.bitList.count(1)])
    for para in paraList:
      paraNum1s.append([iteration,para.bitList.count(1)])

    relFitness.append([iteration,hostWins/(COMPETITION_SIZE*POP_SIZE)])

    #normalise scores and calculate fitness of parasites
    maxScore = max([para.score for para in paraList])
    if maxScore > 0:
      for para in paraList:
        para.score = float(para.score)/maxScore
        para.fitness = ((2.0 * para.score) / (virulence)) - ((para.score * para.score) / (virulence * virulence))

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

def initLists(method=1):
  if method == 1 :
    #initialise lists of participants based on their biases
    hostList = []
    paraList = []
    for x in range(POP_SIZE):
      hostList += [Participant(HOST_BIAS, HOST_MUTATE_CHANCE, [])]
      paraList += [Participant(PARA_BIAS, PARA_MUTATE_CHANCE, [])]
    return (hostList, paraList)
  elif method == 2:
    #initialise lists of participants, parasites starting from all 0s
    hostList = []
    paraList = []
    for x in range(POP_SIZE):
      hostList += [Participant(HOST_BIAS, HOST_MUTATE_CHANCE, [])]
      paraList += [Participant(PARA_BIAS, PARA_MUTATE_CHANCE, [0]*BIT_LENGTH)]
    return (hostList, paraList)


if USE_SEED is None:
  USE_SEED = randrange(sys.maxsize)
  rng = Random(USE_SEED)
print("Seed was:", USE_SEED)

#the different virulence values being used
virulences = [0.5, 0.75, 1.0]

#define figure's size
fig = plt.figure(figsize=[10,12])
#split figure into len(vir) columns
gs0 = fig.add_gridspec( 1, len(virulences), wspace=0.05)

#fill each column
for x in range(len(virulences)):

  #reset all collection lists and seed
  seed(USE_SEED)  
  hostNum1s           = [] #[[iteration, num1s]]
  paraNum1s           = [] #[[iteration, num1s]]
  hostAbsHiffScores   = [] #[[iteration, HiffScore]]                    # actual hiff scores 
  paraAbsHiffScores   = [] #[[iteration, HiffScore]]                    # actual hiff scores 
  hostRelHiffScores   = [] #[[iteration, HiffScore/maxScore]]           # percentage of possible hiff score for given mask
  paraRelHiffScores   = [] #[[iteration, HiffScore/maxScore]]           # percentage of possible hiff score for given mask
  connectednesses     = [] #[[iteration, connectedness]]
  relFitness          = [] #[[iteration, relFitness]]                   # percentage of competitions won by host

  #initialize lists of participants
  (hostList, paraList) = initLists(INIT_METHOD)
  #do iterations with current virulence
  (hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, virulences[x])

  #divide up gridspec and assign array of sublplots
  gs = gridspec.GridSpecFromSubplotSpec(5,1, subplot_spec=gs0[x], hspace=0.05)
  ax = []
  ax.append(plt.subplot(gs[0], xlim=[0,GENERATIONS], xticks=[], ylim=[0,1]))
  ax.append(plt.subplot(gs[1], xlim=[0,GENERATIONS], xticks=[], ylim=[0,maxFitness([1]*BIT_LENGTH)]))
  ax.append(plt.subplot(gs[2], xlim=[0,GENERATIONS], xticks=[], ylim=[0,1]))
  ax.append(plt.subplot(gs[3], xlim=[0,GENERATIONS], xticks=[], ylim=[0,maxFitness([1]*BIT_LENGTH)]))
  ax.append(plt.subplot(gs[4], xlim=[0,GENERATIONS], ylim=[0,1]))

  #plot num1s
  ax[0].title.set_text("Virulence = " + str(virulences[x]))
  if x == 0: 
    ax[0].set_ylabel("Proportion of \n1s in genotype".title())
  else:
    ax[0].yaxis.set_visible(False)
  ax[0].plot([a[0] for a in hostNum1s], [a[1]/BIT_LENGTH for a in hostNum1s], 'o', color='red', markersize=0.1);
  ax[0].plot([a[0] for a in paraNum1s], [a[1]/BIT_LENGTH for a in paraNum1s], 'o', color='blue', markersize=0.1);

  #plot abs hiff
  if x == 0: 
    ax[1].set_ylabel("Absolute HIFF Scores".title())
  else:
    ax[1].yaxis.set_visible(False)
  ax[1].plot([a[0] for a in hostAbsHiffScores], [a[1] for a in hostAbsHiffScores], 'o', color='red', markersize=0.1);
  ax[1].plot([a[0] for a in paraAbsHiffScores], [a[1] for a in paraAbsHiffScores], 'o', color='blue', markersize=0.1);

  #plot relative hiff
  if x == 0: 
    ax[2].set_ylabel("Proportion of possible \nHIFF score".title())
  else:
    ax[2].yaxis.set_visible(False)
  ax[2].plot([a[0] for a in hostRelHiffScores], [a[1] for a in hostRelHiffScores], 'o', color='red', markersize=0.1, clip_on=False, zorder=100);
  ax[2].plot([a[0] for a in paraRelHiffScores], [a[1] for a in paraRelHiffScores], 'o', color='blue', markersize=0.1, clip_on=False, zorder=100);

  #plot conectedness
  if x == 0: 
    ax[3].set_ylabel("mask connectedness".title())
  else:
    ax[3].yaxis.set_visible(False)
  ax[3].plot([a[0] for a in connectednesses], [a[1] for a in connectednesses], 'o', color='red', markersize=0.1, clip_on=False, zorder=100);

  #plot relative fitness
  if x == 0: 
    ax[4].set_ylabel("proportion of \nhost victories".title())
  else:
    ax[4].yaxis.set_visible(False)
  ax[4].plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1)

#display graph
plt.show()
