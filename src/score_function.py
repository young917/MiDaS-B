from collections import defaultdict
import argparse
import pandas as pd
import math
import numpy as np
import sys
import os
import random
from sklearn.linear_model import LinearRegression
from sklearn.metrics import explained_variance_score, max_error, mean_absolute_error, mean_squared_error, mean_squared_log_error, median_absolute_error, r2_score

# DIFF
MAXDIFF = 10000

def read_dist(algoname, dataname, evalname, repeat="2", accuracy=500):
    names = ["0", "1"]
    if (evalname == "sv") or (evalname == "effdiameter"):
        names += ["2"]
    if algoname == "answer":
        path = '../results/answer/' + dataname + "/" + evalname + '.txt'
    else:
        path = '../results/' + algoname + '/' + dataname + "/" + str(repeat) + "/" + evalname + '.txt'
    if os.path.isfile(path) is False:
        return [0]
    
    d = pd.read_csv(path, header=None, index_col=False, names=names)
    if evalname == "intersection":
        ret_dist = list(d["0"] / d["1"])
    elif (evalname == "sv") or (evalname == "effdiameter"):
        ret_dist = list(d["2"])
    elif evalname == 'densification':
        ret_dist = list(d["1"] / d["0"])
    elif evalname == "sizewcc":
        ret_dist = list(d["1"] / max(d["1"]))
    else:
        ret_dist = list(d["1"])

    return ret_dist

def fenalty_score(algoname, dataname, repeat="2", accuracy=500):
    # ESRandom
    compare_repeat_index = repeat
    esrandom_dsf = read_dist("es/add_global_deg_min_0.0000", dataname, "densification", compare_repeat_index, accuracy)
    esrandom_wcc = read_dist("es/add_global_deg_min_0.0000", dataname, "sizewcc", compare_repeat_index, accuracy)

    # Algorithm
    dsf = read_dist(algoname, dataname, "densification", repeat, accuracy)
    wcc = read_dist(algoname, dataname, "sizewcc", repeat, accuracy)
    
    assert len(dsf) == len(esrandom_dsf), str(len(dsf)) + "  <->  " + str(len(esrandom_dsf))
    assert len(wcc) == len(esrandom_wcc), str(len(wcc)) + "  <->  " + str(len(esrandom_wcc))
    assert len(dsf) == len(wcc)

    dsf_diff = 0
    wcc_diff = 0
    for i in range(len(dsf)):
        dsf_diff += (esrandom_dsf[i] - dsf[i])   
        wcc_diff += (wcc[i] - esrandom_wcc[i])
    dsf_diff /= len(dsf)
    wcc_diff /= len(wcc)
    if dsf_diff > 0:
        return True
    elif wcc_diff > 0:
        return True
    else:
        return False

def linscore(algoname, dataname, repeat="2", accuracy=500):
    usecoldict = {
        "intersection" : [0,1],
        "size_avg" : [0,1]
    }
    total_score = 0
    outputdir = '../results/' + algoname + "/" + dataname + "/" 
    if "answer" not in algoname:
        outputdir += repeat + "/"

    for propertyname in ["intersection", "size_avg"]:
        if propertyname in ["intersection"]:
            d = pd.read_csv(outputdir + propertyname + ".txt", names=["Y", "X"], usecols=usecoldict[propertyname])
        else:
            d = pd.read_csv(outputdir + propertyname + ".txt", names=["X", "Y"], usecols=usecoldict[propertyname])
        X = np.array(d["X"])
        Y = np.array(d["Y"])
        X = np.log2(X + 1).reshape(-1,1)
        Y = np.log2(Y + 1).reshape(-1,1)
        reg = LinearRegression().fit(X, Y)
        predY = reg.predict(X)
        total_score += mean_absolute_error(Y, predY)

    return total_score


parser = argparse.ArgumentParser()
parser.add_argument('--dataset', required=False, default="email-Enron-full")
parser.add_argument('--datasetdir', required=False, default="../../dataset/")
parser.add_argument('--algorithm', required=False, default='es/global_deg_min_0.0000')
parser.add_argument('--repeat', required=False, type=str, default="2")
parser.add_argument('--accuracy', required=False, type=int, default=500)
args = parser.parse_args()

outputdir = '../results/' + args.algorithm + "/" + args.dataset + "/" 
if "answer" not in args.algorithm:
    outputdir += args.repeat + "/"
outputfname = outputdir + "score.txt"

diffest = 0
if fenalty_score(args.algorithm, args.dataset, repeat=args.repeat, accuracy=args.accuracy):
    diffest = MAXDIFF
else:
    diffest += linscore(args.algorithm, args.dataset, repeat=args.repeat, accuracy=args.accuracy)

with open(outputfname, "w") as f:
	f.write(str(diffest) + "\n")
    