from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import pylab

BIT_LENGTH       = 128
POP_SIZE         = 50
HOST_BIAS        = 0.5
PARA_BIAS        = 0.5
MUTATE_CHANCE    = 0.03
GENERATIONS      = 600
USE_SEED         = None
VIRULENCE        = 1
RESOURCE_SHARING = True


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
seed(USE_SEED)

someList = [1] * 64

print(recursiveFitness(someList))