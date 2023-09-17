import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import stats
import matplotlib.pyplot as plt
import os
import math
import argparse
plt.rcParams.update({'font.size': 20})

evallist = ["clusteringcoef", "densification", "effdiameter", "overlapness"]
evaldistlist = ["degree", "intersection", "pairdeg", "size", "sv", "wcc"]

ss1 = [0.0, 0.25 ,0.5, 1.0, 2.0]
ss2 = [-1.0, -0.5, -0.25, 0.0, 0.25 ,0.5, 1.0]

searching = {
    "midas_oracle": ["es/add_global_deg_min_%.4f" % (alpha) for alpha in ss1],
    "midasB_oracle": ["essz/add_global_deg_min_%.2f_%.2f" % (alpha, beta) for alpha in ss1 for beta in ss2],
    "baseline": [
                "ns/add_global_deg_0.0000", # RNS
                "ns/add_global_deg_1.0000", # RDN
                "ff/ff_c_0.51_0.20", # FF
                "rw/rw_c_1.00_0.30", # RW
                "ns/add_ordered", # ONS
                "es/add_global_deg_min_0.0000", # RHS
                "tihs", # TIHS
                "midas_oracle" # MiDaS w. Oracle
                 ],
    "final": [
                "ns/add_global_deg_0.0000", # RNS
                "ns/add_global_deg_1.0000", # RDN
                "rw/rw_c_1.00_0.30", # RW
                "ff/ff_c_0.51_0.20",  # FF
                "ns/add_ordered", # ONS
                "essz/add_global_deg_min_0.00_0.00", # RHS
                "tihs", # TIHS
                "midas", # MiDaS w/o. Oracle
                "midasB", # MiDaS-B
    ],
    "ablation": [
                "midas_oracle", # MiDaS w. Oracle
                "midasB_oracle", # MiDaS-B w. Oracle
                "essz/add_global_deg_min_0.25_0.25", # top 1
                "essz/add_global_deg_min_0.00_0.00", # top 5
                "essz/add_global_deg_min_0.00_-0.25", # top 10
                "essz/add_global_deg_min_2.00_-1.00", # last
                "midasB"
    ]
}

dataset = ["email-Enron-full", "email-Eu-full",  
           "contact-high-school", "contact-primary-school",
          "NDC-classes-full", "NDC-substances-full", "tags-ask-ubuntu", "tags-math-sx", 
           "threads-ask-ubuntu", "coauth-MAG-History-full", "coauth-MAG-Geology-full"]
dataset = ["email-Eu-full", "contact-primary-school",
            "NDC-substances-full", "tags-ask-ubuntu"]

# --------------------------------------------------------------------------------------------------------
def get_dist(dirpath, portion):
    dist = {}
    
    # degree, intersection, pairdeg, size: k "," v
    for evalname in ["degree", "intersection", "pairdeg", "size"]:
        dist[evalname] = {}
        tmpdist = {}
        assert os.path.isfile(dirpath + evalname + "%.1f.txt" % (portion)), dirpath + evalname + "%.1f.txt" % (portion)
        with open(dirpath + evalname + "%.1f.txt" % (portion), "r") as f:
            for line in f.readlines():
                k, v = line.rstrip().split(",")
                k, v = float(k), float(v)
                if k > 0:
                    tmpdist[k] = v
        sortedkeys = sorted(list(tmpdist.keys()))
        total = 0
        for k in sortedkeys:
            total += tmpdist[k]
        accum = 0
        for k in sortedkeys:
            accum += (tmpdist[k] / total)
            dist[evalname][k] = accum
            if accum > 1.0000001:
                print(accum)
                print(tmpdist)
            assert accum <= 1.0000001, evalname
    
    # svdist_: k " " v
    dist["sv"] = {}
    assert os.path.isfile(dirpath + "svdist_%.1f.txt" % (portion)), dirpath + "svdist_%.1f.txt" % (portion)
    if os.path.isfile(dirpath + "svdist_%.1f.txt" % (portion)):
        with open(dirpath + "svdist_%.1f.txt" % (portion), "r") as f:
            for line in f.readlines():
                k, v = line.rstrip().split(" ")
                k, v =  float(k), float(v)
                assert v <= 1.0000001, "sv"
                dist["sv"][k] = v

    # wcc: size/numNodes
    dist["wcc"] = {}
    assert os.path.isfile(dirpath + "wcc%.1f.txt" % (portion)), dirpath + "wcc%.1f.txt" % (portion)
    with open(dirpath + "wcc%.1f.txt" % (portion), "r") as f:
        tmplist = []
        for line in f.readlines():
            val = float(line.rstrip())
            tmplist.append(val)
        tmplist = sorted(tmplist, reverse=True)
        total = sum(tmplist)
        accum = 0
        for vi, v in enumerate(tmplist):
            accum += (v / total)
            dist["wcc"][vi] = accum
            if accum > 1.0000001:
                print(accum)
                print(tmpdist)
            assert accum <= 1.0000001, "wcc"
            
    return dist

def get_value(dirpath, portion):
    dist = {}
    
    # clusteringcoef: numhedges - cc  - numnodes
    # effdiameter: numnodes - numhedges - diam
    # overlapness: numhedges - overlapness
    colindex = {
        "clusteringcoef" : 1,
        "effdiameter" : 2,
        "overlapness" : 1
    }
    for evalname in ["clusteringcoef", "effdiameter", "overlapness"]:
        col = colindex[evalname]
        tmplist = []
        totalline = 0
        assert os.path.isfile(dirpath + evalname + ".txt"), dirpath + evalname + ".txt"
        with open(dirpath + evalname + ".txt", "r") as f:
            for line in f.readlines():
                v = float(line.rstrip().split(",")[col])
                tmplist.append(v)
                totalline += 1
        index  = int(portion / (1.0/totalline))
        dist[evalname] = tmplist[index]
    
    # densification ~ density: numnodes - numhedges
    assert os.path.isfile(dirpath + "densification.txt")
    with open(dirpath + "densification.txt" , "r") as f:
        tmplist = []
        totalline = 0
        for line in f.readlines():
            nodes, hedges = line.rstrip().split(",")
            nodes, hedges = int(nodes), int(hedges)
            tmplist.append(hedges / nodes)
            totalline += 1
        index  = int(portion / (1.0/totalline))
        dist["densification"] = tmplist[index]
    
    return dist

def get_dstat(data_dict, answer_dict):
    keyset = set(list(data_dict.keys()) + list(answer_dict.keys()))
    keys = [k for k in keyset]
    keys = sorted(keys)

    answer_cumulsum = 0.0
    data_cumulsum = 0.0
    stat = 0.0
    for k in keys:
        if k in answer_dict:
            answer_cumulsum = answer_dict[k]
            if answer_cumulsum > 1.0001:
                print(answer_dict)
            assert answer_cumulsum <= 1.0001, str(answer_cumulsum)
        if k in data_dict:
            data_cumulsum = data_dict[k]
            assert data_cumulsum <= 1.0001, str(data_cumulsum)
        if stat < math.fabs(answer_cumulsum - data_cumulsum):
            stat = math.fabs(answer_cumulsum - data_cumulsum)
    
    return stat


# --------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--search_name', required=True, type=str, default="baseline")
parser.add_argument('--repeat', type=str, default="2")
args = parser.parse_args()

args.repeat = args.repeat.split(",")
search_name = args.search_name
algonames = searching[search_name]

print(algonames)
print("Read properties")
for portion in [0.1, 0.3, 0.5, 0.7, 0.9]:
    for dataname in dataset:
        evaluation = {}
        alldist = {}
        
        # answer
        alldist["answer"] = {}
        dirpath = '../results/answer/{}/'.format(dataname)
        ret = get_dist(dirpath, portion)
        for e in ret.keys():
            alldist["answer"][e] = ret[e]
        ret = get_value(dirpath, portion)
        for e in ret.keys():
            alldist["answer"][e] = ret[e]
        
        for algoname in algonames:
            # alldist[algoname] = defaultdict(list)
            evaluation[algoname] = defaultdict(list)
            for repeat in args.repeat:
                dirpath = '../results/{}/{}/{}/'.format(algoname, dataname, repeat)
                ret = get_dist(dirpath, portion)
                for e in ret.keys():
                    # alldist[algoname][e].append(ret[e])
                    dstat = get_dstat(ret[e], alldist["answer"][e])
                    evaluation[algoname][e].append(dstat)
                    
                ret = get_value(dirpath, portion)
                for e in ret.keys():
                    # alldist[algoname][e].append(ret[e])
                    diff = abs(ret[e] - alldist["answer"][e])
                    evaluation[algoname][e].append(diff)
            for e in evallist + evaldistlist:
                evaluation[algoname][e] = np.mean(evaluation[algoname][e])
                
        if os.path.isdir("../analyze/csvs/{}/".format(dataname)) is False:
            os.makedirs("../analyze/csvs/{}/".format(dataname))
        with open("../analyze/csvs/%s/eval_%.1f_%s.txt" % (dataname, portion, search_name), "w") as f:
            f.write("algoname")
            for evalname in evallist + evaldistlist:
                f.write("," + evalname)
            f.write("\n")
            for algoname in algonames:
                f.write(algoname)
                for evalname in evallist + evaldistlist:
                    f.write("," + str(evaluation[algoname][evalname]))
                f.write("\n")

# Z-Score & Ranking --------------------------------------------------------------------------
print("Z-Score and Ranking")
for portion in [0.1, 0.3, 0.5, 0.7, 0.9]:     
    for dataname in dataset:   
        d = pd.read_csv("../analyze/csvs/%s/eval_%.1f_%s.txt" % (dataname, portion, search_name))
        for col in (evallist + evaldistlist):
            if d[col].std() != 0:
                d[col] = (d[col] - d[col].mean()) / d[col].std()
        selected_columns = evallist + evaldistlist
        d['avg'] = d[selected_columns].mean(axis=1)
        d.to_csv("../analyze/csvs/%s/eval_diff_norm_%.1f_%s.txt" % (dataname, portion, search_name), index=False)

    for dataname in dataset:
        d = pd.read_csv("../analyze/csvs/%s/eval_%.1f_%s.txt" % (dataname, portion, search_name))
        for ename in (evallist + evaldistlist):
            d[ename] = d[ename].abs().rank(method='min')
        selected_columns = evallist + evaldistlist
        ranks = d[selected_columns]
        d['avg'] = ranks.mean(axis=1)
        d.to_csv("../analyze/csvs/%s/eval_diff_rank_%.1f_%s.txt" % (dataname, portion, search_name), index=False)
    
    for dataname in dataset:
        d = pd.read_csv("../analyze/csvs/%s/eval_%.1f_%s.txt" % (dataname, portion, search_name))
        selected_columns = evallist + evaldistlist
        d['avg'] = d[selected_columns].mean(axis=1)
        d.to_csv("../analyze/csvs/%s/eval_diff_%.1f_%s.txt" % (dataname, portion, search_name), index=False)

# Average over Datasets -------------------------------------------------------------------------
print("Average over Datasets")
for portion in [0.1, 0.3, 0.5, 0.7, 0.9]: 
    result_diff_rank = defaultdict(dict)
    result_diff_norm = defaultdict(dict)
    result_diff = defaultdict(dict)
    for dataname in dataset:
        # diff - rank
        d = pd.read_csv("../analyze/csvs/%s/eval_diff_rank_%.1f_%s.txt" % (dataname, portion, search_name))
        for i, row in d.iterrows():
            algoname = row["algoname"]
            for evalname in evallist + evaldistlist + ["avg"]:
                if evalname not in result_diff_rank[algoname]:
                    result_diff_rank[algoname][evalname] = []
                result_diff_rank[algoname][evalname].append(row[evalname])
        
        # diff - norm
        d = pd.read_csv("../analyze/csvs/%s/eval_diff_norm_%.1f_%s.txt" % (dataname, portion, search_name))
        for i, row in d.iterrows():
            algoname = row["algoname"]
            for evalname in evallist + evaldistlist + ["avg"]:
                if evalname not in result_diff_norm[algoname]:
                    result_diff_norm[algoname][evalname] = []
                result_diff_norm[algoname][evalname].append(row[evalname])
        
        # diff
        d = pd.read_csv("../analyze/csvs/%s/eval_diff_%.1f_%s.txt" % (dataname, portion, search_name))
        for i, row in d.iterrows():
            algoname = row["algoname"]
            for evalname in evallist + evaldistlist + ["avg"]:
                if evalname not in result_diff[algoname]:
                    result_diff[algoname][evalname] = []
                result_diff[algoname][evalname].append(row[evalname])

    for algoname in result_diff_rank.keys():
        for evalname in evallist + evaldistlist + ["avg"]:
            assert len(result_diff_rank[algoname][evalname]) == len(dataset), len(result_diff_rank[algoname][evalname])
            if evalname == "avg":
                result_diff_rank[algoname][evalname + " std"] = np.std(np.array(result_diff_rank[algoname][evalname]))
            result_diff_rank[algoname][evalname] = np.mean(np.array(result_diff_rank[algoname][evalname]))
            
            assert len(result_diff_norm[algoname][evalname]) == len(dataset), len(result_diff_norm[algoname][evalname])
            if evalname == "avg":
                result_diff_norm[algoname][evalname + " std"] = np.std(np.array(result_diff_norm[algoname][evalname]))
            result_diff_norm[algoname][evalname] = np.mean(np.array(result_diff_norm[algoname][evalname]))
            
            assert len(result_diff[algoname][evalname]) == len(dataset), len(result_diff[algoname][evalname])
            if evalname == "avg":
                result_diff[algoname][evalname + " std"] = np.std(np.array(result_diff[algoname][evalname]))
            result_diff[algoname][evalname] = np.mean(np.array(result_diff[algoname][evalname]))

    with open("../analyze/csvs/eval_diff_rank_all_%.1f_%s.txt" % (portion,search_name), "w") as f:
        f.write(",".join(["algoname"] + evallist + evaldistlist + ["avg", "avg std"]) + "\n")
        for algoname in result_diff_rank.keys():
            f.write(algoname)
            for evalname in evallist + evaldistlist + ["avg", "avg std"]:
                f.write("," + str(result_diff_rank[algoname][evalname]))
            f.write("\n")

    with open("../analyze/csvs/eval_diff_norm_all_%.1f_%s.txt" % (portion,search_name), "w") as f:
        f.write(",".join(["algoname"] + evallist + evaldistlist + ["avg", "avg std"]) + "\n")
        for algoname in result_diff_norm.keys():
            f.write(algoname)
            for evalname in evallist + evaldistlist + ["avg", "avg std"]:
                f.write("," + str(result_diff_norm[algoname][evalname]))
            f.write("\n")
    
    with open("../analyze/csvs/eval_diff_all_%.1f_%s.txt" % (portion,search_name), "w") as f:
        f.write(",".join(["algoname"] + evallist + evaldistlist + ["avg", "avg std"]) + "\n")
        for algoname in result_diff.keys():
            f.write(algoname)
            for evalname in evallist + evaldistlist + ["avg", "avg std"]:
                f.write("," + str(result_diff[algoname][evalname]))
            f.write("\n")


# Aggregate Portion ---------------------------------------------------------------------------
result_diff_rank = defaultdict(dict)
result_diff_norm = defaultdict(dict)
result_diff = defaultdict(dict)
for portion in [0.1, 0.3, 0.5, 0.7, 0.9]:
    # diff - rank
    d = pd.read_csv("../analyze/csvs/eval_diff_rank_all_%.1f_%s.txt" % (portion, search_name))
    for i, row in d.iterrows():
        algoname = row["algoname"]
        for evalname in  evaldistlist + evallist + ["avg"]:
            if evalname not in result_diff_rank[algoname]:
                result_diff_rank[algoname][evalname] = []
            result_diff_rank[algoname][evalname].append(row[evalname])
    
    # diff - norm
    d = pd.read_csv("../analyze/csvs/eval_diff_norm_all_%.1f_%s.txt" % (portion, search_name))
    for i, row in d.iterrows():
        algoname = row["algoname"]
        for evalname in evaldistlist + evallist + ["avg"]:
            if evalname not in result_diff_norm[algoname]:
                result_diff_norm[algoname][evalname] = []
            result_diff_norm[algoname][evalname].append(row[evalname])

    # diff
    d = pd.read_csv("../analyze/csvs/eval_diff_all_%.1f_%s.txt" % (portion, search_name))
    for i, row in d.iterrows():
        algoname = row["algoname"]
        for evalname in evaldistlist + evallist + ["avg"]:
            if evalname not in result_diff[algoname]:
                result_diff[algoname][evalname] = []
            result_diff[algoname][evalname].append(row[evalname]) 

result_diff_rank_all = {}
result_diff_norm_all = {}
result_diff_all = {}

for algoname in result_diff_rank.keys():
    result_diff_rank_all[algoname] = {}
    result_diff_norm_all[algoname] = {}
    result_diff_all[algoname] = {}
    for evalname in evaldistlist + evallist + ["avg"]:
        result_diff_rank_all[algoname][evalname] = np.mean(np.array(result_diff_rank[algoname][evalname]))
        result_diff_norm_all[algoname][evalname] = np.mean(np.array(result_diff_norm[algoname][evalname]))
        result_diff_all[algoname][evalname] = np.mean(np.array(result_diff[algoname][evalname]))

with open("../analyze/csvs/eval_diff_rank_all_%s.txt" % (search_name), "w") as f:
    f.write(",".join(["algoname"] + evaldistlist + evallist + ["avg"]) + "\n")
    for algoname in searching[search_name]: #result_diff_rank_all.keys():
        f.write(algoname)
        for evalname in evaldistlist + evallist + ["avg"]:
            f.write("," + str(result_diff_rank_all[algoname][evalname]))
        f.write("\n")

with open("../analyze/csvs/eval_diff_norm_all_%s.txt" % (search_name), "w") as f:
    f.write(",".join(["algoname"] + evaldistlist + evallist + ["avg"]) + "\n")
    for algoname in searching[search_name]: #result_diff_norm_all.keys():
        f.write(algoname)
        for evalname in evaldistlist + evallist + ["avg"]:
            f.write("," + str(result_diff_norm_all[algoname][evalname]))
        f.write("\n")

with open("../analyze/csvs/eval_diff_all_%s.txt" % (search_name), "w") as f:
    f.write(",".join(["algoname"] + evaldistlist + evallist + ["avg"]) + "\n")
    for algoname in searching[search_name]: #result_diff_norm_all.keys():
        f.write(algoname)
        for evalname in evaldistlist + evallist + ["avg"]:
            f.write("," + str(result_diff_all[algoname][evalname]))
        f.write("\n")

print("End ", search_name)
