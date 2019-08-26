#
# H-IFF simple evolutionary algorithm
#

from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import pylab
import statistics
import multiprocessing

BIT_LENGTH       = 32
POP_SIZE         = 25
HOST_BIAS        = 0.5
MUTATE_CHANCE    = 0.03
GENERATIONS      = 600
USE_SEED         = None
SELECTION_SIZE   = 5


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


def mainLoop(nIters, iterOffset):
  global hostAbsHiffScores, hostRelHiffScores, hostNum1s
  global paraAbsHiffScores, paraRelHiffScores, paraNum1s

  hostList = []
  for x in range(POP_SIZE):
    hostList   += [Participant(HOST_BIAS)]

  for iteration in range(iterOffset, nIters):
    #mutate all participants and recalculate scores
    for host in hostList:
      host.mutate()
      hostNum1s.append([iteration,host.bitList.count(1)])
      host.fitness = recursiveFitness(host.bitList)
      hostAbsHiffScores.append([iteration,host.fitness])

    #asexual breeeding with tournement size 5
    newHostList = []
    for index in range(POP_SIZE):
      hostSample = [hostList[i] for i in choices(range(POP_SIZE),k=SELECTION_SIZE)] #use choices for replacement
      best = hostSample[0]
      for host in hostSample:
        best = host if host.fitness > best.fitness else best
      newHostList.append(copy.deepcopy(best))
    hostList = newHostList

  print(max([recursiveFitness(host.bitList) for host in hostList]))
  return hostList

def initList():  
  #initialise list of participants
  hostList = []
  for x in range(POP_SIZE):
    hostList   += [Participant(HOST_BIAS)]
  return hostList 


if USE_SEED is None:
  USE_SEED = randrange(sys.maxsize)
  rng = Random(USE_SEED)
print("Seed was:", USE_SEED)
seed(USE_SEED)

hostList = []

hostNum1s           = [] #[[iteration, num1s]]
hostAbsHiffScores   = [] #[[iteration, HiffScore]]                    # actual hiff scores 

##############use this for running many times###############
if __name__ == '__main__':
  for n in range(20):
    jobs = []
    for x in range(50):
      p = multiprocessing.Process(target=mainLoop, args=(GENERATIONS,0))
      jobs.append(p)
      p.start()
    p.join()

##############use this for running once and graphing###############

# #do 600 iterations
# hostList = mainLoop(GENERATIONS, 0)

# #define sublot parameters
# f, (ax0, ax1) = plt.subplots(2, 1, sharex=True, subplot_kw={'xlim': {0,GENERATIONS}})

# #plot num1s
# ax0.title.set_text("Proportion of 1s in genotype")
# ax0.set_ylim([0,1])
# ax0.plot([a[0] for a in hostNum1s], [a[1]/BIT_LENGTH for a in hostNum1s], 'o', color='red', markersize=0.1);

# #plot abs hiff
# ax1.title.set_text("Absolute HIFF Scores")
# ax1.set_ylim([0,maxFitness([1]*BIT_LENGTH)])
# ax1.plot([a[0] for a in hostAbsHiffScores], [a[1] for a in hostAbsHiffScores], 'o', color='red', markersize=0.1);

# #display graph
# plt.show()
