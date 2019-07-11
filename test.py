from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import pylab

class Participant():
  """Host or Parasite parent"""
  bitList = [False] * 100
  score = 0
  fitness = 0
  bias = 0.5
  victors = []

  def __init__(self, bias):
    self.bitList = []
    for x in range(100):
      self.bitList += [False] if random()<bias else [True]
    self.bias = bias

part = Participant(0.5)

part1 = Participant(0.75)

part.victors += [part1]

part.victors[0].score += 1

print(part1.score)