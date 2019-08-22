from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import pylab
import itertools

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

# resultsList = [] # [[n1s, score]]
# bestScores  = []
# worstScores = []

# for n1s in range(32 + 1):
#   n0s = 32-n1s
#   testList = [0] * n0s + [1] * n1s
#   bestScore = 0
#   worstScore = recursiveFitness([0]*32)
#   for iteration in range(1000):
#     shuffle(testList)
#     score = recursiveFitness(testList)
#     resultsList.append([n1s, score])
#     if score > bestScore:
#       bestScore = score
#     if score < worstScore:
#       worstScore = score
#   bestScores.append([n1s,bestScore])
#   worstScores.append([n1s.worstScore])

# plt.plot([a[0] for a in resultsList], [a[1] for a in resultsList], 'o', color='red', markersize=0.1)
# plt.xlim([0, 32])
# plt.ylim([0, recursiveFitness([0]*32)])

# plt.show()
# for item in list(itertools.product([0, 1], repeat=4)):
#   print(item)
#   print(recursiveFitness(item))

def mutate(bitList):
  newBitList = []
  for bit in bitList:
    if random() < 0.03:
      newBitList += [1] if random()<0.5 else [0]
    else:
      newBitList += [bit]
  return newBitList

lst = [0,0,1,1]
trials = 10000000
counter = 0

for x in range(trials):
  if mutate(lst) == [0,0,0,0] or mutate(lst) == [1,1,1,1]:
    counter+=1

print(trials/counter)

