import os
from collections import defaultdict
import numpy as np
import pandas as pd

dataset = [
    "email-Enron-full", "email-Eu-full",
    "contact-high-school", "contact-primary-school",
    "NDC-classes-full", "NDC-substances-full",
    "tags-ask-ubuntu", "tags-math-sx",
    "threads-ask-ubuntu", "coauth-MAG-History-full", "coauth-MAG-Geology-full"
]
dataset = [
    "email-Eu-full", "contact-primary-school", "NDC-substances-full",
    "tags-ask-ubuntu", "coauth-MAG-History-full"]
ASS = ["0.00", "0.25", "0.50", "1.00", "2.00"]
BSS = ["-1.00", "-0.50", "-0.25", "0.00", "0.25", "0.50", "1.00"]

# repeat = [1,2,3]
repeat = [1]

data2bestalgo = defaultdict()

for dataname in dataset:
    score_mat = np.zeros((len(ASS), len(BSS)))
    for aindex, alpha in enumerate(ASS):
        for bindex, beta in enumerate(BSS):
            reject_times = 0
            for repeatindex in repeat:
                fname = "../results/essz/add_global_deg_min_{}_{}/{}/{}/score.txt".format(alpha, beta, dataname, repeatindex)
                with open(fname, "r") as f:
                    score = float(f.readline().rstrip())
                score_mat[aindex][bindex] += score
            score_mat[aindex][bindex] /= len(repeat)
    
    ind0, ind1 = np.unravel_index(np.argmin(score_mat, axis=None), score_mat.shape)
    data2bestalgo[dataname] = ASS[ind0] + "_" + BSS[ind1]

name = "midasB"
dirpath = "../results/" + name + "/"

if os.path.isdir(dirpath) is False:
    os.makedirs(dirpath)
if os.path.isfile(dirpath + "search_result.txt"):
    os.remove(dirpath + "search_result.txt")

with open(dirpath + "search_result.txt", "w") as f:
    for dataname in dataset:
        param = data2bestalgo[dataname]    
        line = "/home/dmlab/minyoung/BackInTime_Sampling_cpp/results/essz/add_global_deg_min_{}/{}/".format(param, dataname) + "\n"
        f.write(line)