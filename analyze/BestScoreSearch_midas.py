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
repeatset = range(1,4)
ASS = ["0.0000", "0.2500", "0.5000", "1.0000", "2.0000"]

data2bestalgo = defaultdict()

for dataname in dataset:
    scores = []
    for aindex, alpha in enumerate(ASS):
        curscore = 0
        for repeatindex in repeatset:
            fname = "../results/es/add_global_deg_min_{}/{}/{}/score.txt".format(alpha, dataname, repeatindex)
            with open(fname, "r") as f:
                score = float(f.readline().rstrip())
            curscore += score
        curscore /= len(repeatset)
        scores.append(curscore)
    
    str_score = []
    for score in scores:
        str_score.append("%.4f" % (score))
    print(dataname + "   " + ", ".join(str_score))

    best_aindex = np.argmin(np.array(scores))
    data2bestalgo[dataname] = ASS[best_aindex]

name = "midas"
dirpath = "../results/" + name + "/"

if os.path.isdir(dirpath) is False:
    os.makedirs(dirpath)
if os.path.isfile(dirpath + "search_result.txt"):
    os.remove(dirpath + "search_result.txt")

with open(dirpath + "search_result.txt", "w") as f:
    for dataname in dataset:
        param = data2bestalgo[dataname]
        line = "/home/dmlab/minyoung/BackInTime_Sampling_cpp/results/es/add_global_deg_min_{}/{}/".format(param, dataname) + "\n"
        f.write(line)
