import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import stats
import os
import math
import argparse

evallist = ["clusteringcoef", "densification", "effdiameter", "overlapness"]
evaldistlist = ["degree", "intersection", "pairdeg", "size", "sv", "wcc"]

dataset = ["email-Enron-full", "email-Eu-full",  
           "contact-high-school", "contact-primary-school",
          "NDC-classes-full", "NDC-substances-full", "tags-ask-ubuntu", "tags-math-sx", 
           "threads-ask-ubuntu", "coauth-MAG-History-full", "coauth-MAG-Geology-full"]
dataset = ["email-Eu-full", "contact-primary-school",
            "NDC-substances-full", "tags-ask-ubuntu"]

parser = argparse.ArgumentParser()
parser.add_argument('--search_name', required=True)
args = parser.parse_args()

# search_name = "midas_oracle"
# search_name = "midasB_oracle"

for dataname in dataset:
    result_diff_rank = defaultdict(dict)
    result_diff_norm = defaultdict(dict)
    for portion in [0.1, 0.3, 0.5, 0.7, 0.9]:
        # diff - rank
        d = pd.read_csv("csvs/%s/eval_diff_rank_%.1f_%s.txt" % (dataname, portion, args.search_name))
        for i, row in d.iterrows():
            algoname = row["algoname"]
            for evalname in evallist + evaldistlist + ["avg"]:
                if evalname not in result_diff_rank[algoname]:
                    result_diff_rank[algoname][evalname] = []
                result_diff_rank[algoname][evalname].append(row[evalname])

        # diff - norm
        d = pd.read_csv("csvs/%s/eval_diff_norm_%.1f_%s.txt" % (dataname, portion, args.search_name))
        for i, row in d.iterrows():
            algoname = row["algoname"]
            for evalname in evallist + evaldistlist + ["avg"]:
                if evalname not in result_diff_norm[algoname]:
                    result_diff_norm[algoname][evalname] = []
                result_diff_norm[algoname][evalname].append(row[evalname])
    
    for algoname in result_diff_rank.keys():
        for evalname in evallist + evaldistlist + ["avg"]:
            result_diff_rank[algoname][evalname] = np.mean(np.array(result_diff_rank[algoname][evalname]))
            result_diff_norm[algoname][evalname] = np.mean(np.array(result_diff_norm[algoname][evalname]))
    
    with open("csvs/%s/eval_diff_rank_all_%s.txt" % (dataname, args.search_name), "w") as f:
        f.write(",".join(["algoname"] + evallist + evaldistlist + ["avg"]) + "\n")
        for algoname in result_diff_rank.keys():
            f.write(algoname)
            for evalname in evallist + evaldistlist + ["avg"]:
                f.write("," + str(result_diff_rank[algoname][evalname]))
            f.write("\n")

    with open("csvs/%s/eval_diff_norm_all_%s.txt" % (dataname, args.search_name), "w") as f:
        f.write(",".join(["algoname"] + evallist + evaldistlist + ["avg"]) + "\n")
        for algoname in result_diff_norm.keys():
            f.write(algoname)
            for evalname in evallist + evaldistlist + ["avg"]:
                f.write("," + str(result_diff_norm[algoname][evalname]))
            f.write("\n")

bestalgo = {}
for dataname in dataset:
    result = {}
    
    d = pd.read_csv("csvs/%s/eval_diff_norm_all_%s.txt" % (dataname, args.search_name))
    minnorm = min(d["avg"])
    maxnorm = max(d["avg"])
    for i, row in d.iterrows():
        algoname = row["algoname"]
        norm = (row["avg"] - minnorm) / (maxnorm - minnorm)
        result[algoname] = norm
    
    
    d = pd.read_csv("csvs/%s/eval_diff_rank_all_%s.txt" % (dataname, args.search_name))
    minrank = min(d["avg"])
    maxrank = max(d["avg"])
    for i, row in d.iterrows():
        algoname = row["algoname"]
        rank = (row["avg"] - minrank) / (maxrank - minrank)
        result[algoname] += rank
    
    bestalgoname = ""
    bestalgoval = 1000
    for k, v in result.items():
        if bestalgoval > v:
            bestalgoname = k
            bestalgoval = v
    bestalgo[dataname] = bestalgoname

name = args.search_name
dirpath = "../results/" + name + "/"

if os.path.isdir(dirpath) is False:
    os.makedirs(dirpath)
if os.path.isfile(dirpath + "search_result.txt"):
    os.remove(dirpath + "search_result.txt")

for dataname in dataset:
    with open(dirpath + "search_result.txt", "a") as f:
        algoname = bestalgo[dataname]
        line = "/home/dmlab/minyoung/BackInTime_Sampling_cpp/results/{}/{}/".format(algoname, dataname) + "\n"
        f.write(line)
