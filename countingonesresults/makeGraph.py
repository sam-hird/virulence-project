import sys
import ast
import numpy as np
import matplotlib.pyplot as plt

mutationrates = ["0.005","0.020","0.035","0.050"]
virulences = ["0.5", "1.0", "0.75", "phantom"]
markers = ["s", "o", "x", "+"]

f, ax = plt.subplots(1, 4, sharey=True)
for graph in range(4):
  for dataset in range(len(virulences)):
    #read data from file
    file = open("results_"+ mutationrates[graph] + "_" + virulences[dataset] + ".txt")
    contents = file.read()
    array = ast.literal_eval(contents)

    #add data to graph
    ax[graph].plot([a[0] for a in array] ,[a[1] for a in array], marker=markers[dataset], fillstyle="none", markersize=10, linestyle=":", label="λ="+virulences[dataset] if dataset<3 else virulences[dataset])
    ax[graph].set_xlabel("Parasite Bias")
  ax[graph].legend()
  ax[graph].title.set_text("Mutation rate = " + mutationrates[graph])
ax[0].set_ylabel("Disengaged Generations")

plt.show()
  