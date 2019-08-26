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
POP_SIZE              = 10
HOST_BIAS             = 0.5
PARA_BIAS             = 0.5
HOST_MUTATE_CHANCE    = 0.03
PARA_MUTATE_CHANCE    = 0.03
GENERATIONS           = 1000
USE_SEED              = None
VIRULENCE             = 1
SELECTION_SIZE        = 5
COMPETITION_SIZE      = 5
THRESHOLD             = 0.9
INIT_METHOD           = 1
MUTATE_METHOD         = 3


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
        self.bitList += [1] if random()<0.5 else [0]
    else:
      self.bitList = bitList

    self.mutationRate = mutationRate
    self.bias = bias

  def mutate(self, method):
    if method == 1 :
      newBitList = []
      for bit in self.bitList:
        if random() < self.mutationRate:
          newBitList += [1] if random()<self.bias else [0]
        else:
          newBitList += [bit]
      self.bitList = newBitList

    elif method == 2: 
      newBitList = []
      for bit in range(len(self.bitList)):        
        if random() < self.mutationRate:
          newBitList += ([self.bitList[bit+1]] if bit == 0 else [self.bitList[bit-1]])
        else:
          newBitList += [self.bitList[bit]]
      self.bitList = newBitList

    elif method == 3: 
      newBitList = []
      for bit in range(len(self.bitList)):        
        if random() < self.mutationRate:
          neighbor = (self.bitList[bit+1] if bit == 0 else self.bitList[bit-1])
          newBitList += [neighbor] if random() < self.bias else [1 - neighbor]
        else:
          newBitList += [self.bitList[bit]]
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
  global hostHiffScores, hostProportions, hostNum1s, hostConnectednesses
  global                 paraProportions, paraNum1s, paraConnectednesses
  for iteration in range(iterOffset, nIters):
    #mutate all participants and reset scores
    for host in hostList:
      host.mutate(MUTATE_METHOD)
      host.score = 0
    for para in paraList:
      para.mutate(MUTATE_METHOD)
      para.score = 0

    #shuffle both lists and have the two populations compete pairwise
    hostWins = 0
    for x in range(COMPETITION_SIZE):
      shuffle(hostList)
      for popIndex in range(POP_SIZE):
        #mask lists
        (masked1, masked2) = maskLists(paraList[popIndex].bitList,hostList[popIndex].bitList)
        (mask1 , mask2 ) = maskLists(paraList[popIndex].bitList,[1]*BIT_LENGTH)

        #record mask details
        connectedness1 = 0 if (mask1.count(1) == 0) else recursiveFitness(mask1) / recursiveFitness([1] * mask1.count(1) + [None] * (BIT_LENGTH - mask1.count(1)))
        hostProportion = 0 if recursiveFitness(mask1) == 0 else recursiveFitness(masked1)/recursiveFitness(mask1)

        #lines from mask-complement system
        #connectedness2 = 0 if (mask2.count(1) == 0) else recursiveFitness(mask2) / recursiveFitness([1] * mask2.count(1) + [None] * (BIT_LENGTH - mask2.count(1)))
        #paraProportion = 0 if recursiveFitness(mask2) == 0 else recursiveFitness(masked2)/recursiveFitness(mask2) 

        #calulate winner of this round
        if hostProportion >= THRESHOLD:
          #host wins
          hostList[popIndex].score += 1
          hostWins +=1
        else:
          #para wins
          paraList[popIndex].score += 1

        #store details about this round
        hostConnectednesses.append([iteration,connectedness1])
        hostProportions.append([iteration,hostProportion])

        #lines from mask-complement system
        #paraConnectednesses.append([iteration,connectedness2])
        #paraProportions.append([iteration,paraProportion])


    for host in hostList:
      host.fitness = host.score
      hostHiffScores.append([iteration, recursiveFitness(host.bitList)])
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

    #asexual breeeding
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
fig = plt.figure(figsize=[10,5*len(virulences)])

#split figure into len(vir) columns
gs0 = fig.add_gridspec( 1, len(virulences), wspace=0.05)

#fill each column
for x in range(len(virulences)):

  #reset all collection lists and seed
  seed(USE_SEED)  

  hostHiffScores    = [] #[[iteration, HiffScore]]                    # actual hiff scores 

  hostProportions   = [] #[[iteration, HiffScore/maxScore]]           # percentage of possible hiff score for given mask
  paraProportions   = [] #[[iteration, HiffScore/maxScore]]           # percentage of possible hiff score for given mask

  hostNum1s           = [] #[[iteration, num1s]]
  paraNum1s           = [] #[[iteration, num1s]]

  hostConnectednesses     = [] #[[iteration, connectedness]]
  paraConnectednesses     = [] #[[iteration, connectedness]]

  relFitness          = [] #[[iteration, relFitness]]                   # percentage of competitions won by host

  #initialize lists of participants
  (hostList, paraList) = initLists(INIT_METHOD)
  #do iterations with current virulence
  (hostList, paraList) = mainLoop(GENERATIONS, 0, hostList, paraList, virulences[x])
  print('virulence = ' + str(virulences[x]))
  #show resulting hosts in console
  for host in hostList:
    print(host.bitList);
  print('\n')
  #divide up gridspec and assign array of sublplots
  gs = gridspec.GridSpecFromSubplotSpec(5,1, subplot_spec=gs0[x], hspace=0.05)
  gs.height_ratios = [5,5,5,5,1]
  ax = []
  ax.append(plt.subplot(gs[0], xlim=[0,GENERATIONS], xticks=[], ylim=[0,recursiveFitness([1]*BIT_LENGTH)]))
  ax.append(plt.subplot(gs[1], xlim=[0,GENERATIONS], xticks=[], ylim=[0,1]))
  ax.append(plt.subplot(gs[2], xlim=[0,GENERATIONS], xticks=[], ylim=[0,BIT_LENGTH]))
  ax.append(plt.subplot(gs[3], xlim=[0,GENERATIONS], xticks=[], ylim=[0,1]))
  ax.append(plt.subplot(gs[4], xlim=[0,GENERATIONS], ylim=[0,1]))

  #plot H-IFF scores
  ax[0].title.set_text("Virulence = " + str(virulences[x]))
  if x == 0: 
    ax[0].set_ylabel("H-IFF score of hosts".title())
  else:
    ax[0].yaxis.set_visible(False)
  ax[0].plot([a[0] for a in hostHiffScores], [a[1] for a in hostHiffScores], 'o', color='red', markersize=0.1, clip_on=False, zorder=100);

  #plot proportions
  if x == 0: 
    ax[1].set_ylabel("proportion of possible H-IFF".title())
  else:
    ax[1].yaxis.set_visible(False)
  ax[1].plot([a[0] for a in hostProportions], [a[1] for a in hostProportions], 'o', color='red', markersize=0.1, clip_on=False, zorder=100);
  ax[1].plot([a[0] for a in paraProportions], [a[1] for a in paraProportions], 'o', color='blue', markersize=0.1, clip_on=False, zorder=100);

  #plot num1s
  if x == 0: 
    ax[2].set_ylabel("number of 1s in genotype".title())
  else:
    ax[2].yaxis.set_visible(False)
  ax[2].plot([a[0] for a in hostNum1s], [a[1] for a in hostNum1s], 'o', color='red', markersize=0.1, clip_on=False, zorder=100);
  ax[2].plot([a[0] for a in paraNum1s], [a[1] for a in paraNum1s], 'o', color='blue', markersize=0.1, clip_on=False, zorder=100);

  #plot conectedness
  if x == 0: 
    ax[3].set_ylabel("mask connectedness".title())
  else:
    ax[3].yaxis.set_visible(False)
  ax[3].plot([a[0] for a in hostConnectednesses], [a[1] for a in hostConnectednesses], 'o', color='red', markersize=0.1, clip_on=False, zorder=100);
  ax[3].plot([a[0] for a in paraConnectednesses], [a[1] for a in paraConnectednesses], 'o', color='blue', markersize=0.1, clip_on=False, zorder=100);

  #plot relative fitness
  if x == 0: 
    ax[4].set_ylabel("proportion of \nhost victories".title())
  else:
    ax[4].yaxis.set_visible(False)
  ax[4].set_xlabel("generations".title())
  ax[4].plot([a[0] for a in relFitness], [a[1] for a in relFitness], '.', color="black", markersize=1, clip_on=False, zorder=100)

#display graph
plt.show()
