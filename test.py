from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import pylab

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

def randomList():
  returnList = []
  for x in range(8):
    returnList += 0 if random() < 0.5 else 1
  return returnList

hostLists = [[0,1,0,1,0,1,0,1], [0,0,1,1,0,0,1,1], [0,0,0,0,1,1,1,1], [0,0,0,0,0,0,0,0], [1,1,1,1,1,1,1,1]]



for x in range(20):
  paraLists = []
  for i in range(5):
    paraLists.append((randomList()))
  for listx in hostLists:
    print(recursiveFitness(listx))