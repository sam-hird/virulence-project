from random import *
from operator import attrgetter
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import pylab

someList = [1,2,3,4]

someList[2] = None

print( someList)
print(len(someList))