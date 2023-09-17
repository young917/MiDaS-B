import os
import sys
import argparse
from collections import defaultdict
from tqdm import tqdm
import numpy as np
from itertools import chain
import networkx as nx
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import svds, eigs, norm
from scipy.stats import skew
import snap
import powerlaw
import math
import time

RESULTDIR = "/home/dmlab/minyoung/BackInTime_Sampling_cpp"

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', required=False, default="email-Enron-full")
parser.add_argument('--datasetdir', required=False, default="../../dataset/")
parser.add_argument('--algorithm', required=False, default='es/global_deg_min_0.0000')
parser.add_argument('--repeat', required=False, type=int, default=-1)
parser.add_argument('--accuracy', required=False, type=int, default=500)

parser.add_argument('--degree', required=False, action='store_true')
parser.add_argument('--size', required=False, action='store_true')
parser.add_argument('--pairdegree', required=False, action='store_true')

parser.add_argument('--verbose', required=False, action='store_true')

parser.add_argument('--recalculate', required=False, action='store_true')

args = parser.parse_args()
outputdir = '../results/' + args.algorithm + "/" + args.dataset + "/"
if args.repeat > 0:
    outputdir += str(args.repeat) + "/"
# print(outputdir)
# print(args.effdiam)
# print(args.sv)
print("Start " + outputdir)

datasetdir = "../../dataset/"
entire_nodeset = set([])
entire_number_of_hyperedges = 0
datapath = datasetdir + args.dataset + "_final.txt"
datanode2index = {}
with open(datapath, "r") as f:
    for line in f.readlines():
        hyperedge = line[:-1].split(",")
        entire_number_of_hyperedges += 1
        for n in hyperedge:
            if int(n) not in datanode2index:
                datanode2index[int(n)] = len(datanode2index)
            nodeindex = datanode2index[int(n)]
            entire_nodeset.add(nodeindex)
entire_number_of_nodes = len(entire_nodeset)

hyperedge_steps = []
unit_numhedge = int(math.ceil(entire_number_of_hyperedges / args.accuracy))
for _i in range(entire_number_of_hyperedges):
    if (_i > 0) and (_i % unit_numhedge == 0):
        hyperedge_steps.append(_i)
hyperedge_steps.append(entire_number_of_hyperedges)

# degree ------------------------------------------------------------------------------
if args.degree:
    inputpath = outputdir + "sampled_graph.txt"
    if args.algorithm == "answer":
        inputpath = args.datasetdir + args.dataset + "_final.txt"
    assert(os.path.isfile(inputpath)), inputpath

    runflag = True
    if os.path.isfile(outputdir + "degree_avg.txt"):
        if args.recalculate:
            os.remove(outputdir + "degree_avg.txt")
        else:
            runflag = False
            num_line = 0
            with open(outputdir + "degree_avg.txt", "r") as f:
                for line in f.readlines():
                    h, _ = line.rstrip().split(",")
                    if num_line >= len(hyperedge_steps):
                        runflag = True
                        os.remove(outputdir + "degree_avg.txt")
                        break
                    if int(h) != hyperedge_steps[num_line]:
                        runflag = True
                        os.remove(outputdir + "degree_avg.txt")
                        break
                    num_line += 1
            if num_line < len(hyperedge_steps):
                runflag = True
                os.remove(outputdir + "degree_avg.txt")
            elif num_line == len(hyperedge_steps):
                print("Exist ", outputdir + "degree_avg.txt")
    if runflag:
        step = 0
        node_degree = defaultdict(int)
        number_of_edges = 0
        nodeset = set()

        with open(inputpath, "r") as f:
            for line in f.readlines():
                line = line[:-1]
                hyperedge = [int(_n) for _n in line.split(",")]
                number_of_edges += 1
                for n in hyperedge:
                    node_degree[n] += 1
                    nodeset.add(n)

                if number_of_edges == hyperedge_steps[step]:
                    if args.verbose and step % 10 == 0:
                        print("[" + str(step) + "]")
                    step += 1
                    avgdegree = 0
                    number_of_nodes = len(nodeset)
                    for n in node_degree:
                        deg = node_degree[n]
                        avgdegree += deg / number_of_nodes
                    with open(outputdir + "degree_avg.txt", "+a") as f:
                        f.write(str(number_of_edges) + "," + str(avgdegree) + "\n")
        print("End Degree Avg")

# size ------------------------------------------------------------------------------
if args.size:
    inputpath = outputdir + "sampled_graph.txt"
    if args.algorithm == "answer":
        inputpath = args.datasetdir + args.dataset + "_final.txt"
    assert(os.path.isfile(inputpath)), inputpath

    runflag = True
    if os.path.isfile(outputdir + "size_avg.txt"):
        if args.recalculate:
            os.remove(outputdir + "size_avg.txt")
        else:
            runflag = False
            num_line = 0
            with open(outputdir + "size_avg.txt", "r") as f:
                for line in f.readlines():
                    h, _ = line.rstrip().split(",")
                    if num_line >= len(hyperedge_steps):
                        runflag = True
                        os.remove(outputdir + "size_avg.txt")
                        break
                    if int(h) != hyperedge_steps[num_line]:
                        runflag = True
                        os.remove(outputdir + "size_avg.txt")
                        break
                    num_line += 1
            if num_line < len(hyperedge_steps):
                runflag = True
                os.remove(outputdir + "size_avg.txt")
            elif num_line == len(hyperedge_steps):
                print("Exist ", outputdir + "size_avg.txt")
    if runflag:
        number_of_edges = 0
        step = 0
        sumsize = 0

        start_time = time.process_time()
        with open(inputpath, "r") as f:
            for line in f.readlines():
                line = line[:-1]
                hyperedge = [int(_n) for _n in line.split(",")]
                number_of_edges += 1
                sumsize += len(hyperedge)

                if number_of_edges == hyperedge_steps[step]:
                    if args.verbose and step % 10 == 0:
                        print("[" + str(step) + "]")
                    step += 1
                    avgsize = sumsize / number_of_edges

                    with open(outputdir + "size_avg.txt", "+a") as f:
                        f.write(str(number_of_edges) + "," + str(avgsize) + "\n")
        end_time = time.process_time()
        with open(outputdir + "time_size_avg.txt", "w") as f:
            f.write("{} ms\n".format(int(round((end_time - start_time) * 1000))))
        print("End Size Avg")

# pairdegree ------------------------------------------------------------------------------
if args.pairdegree:
    inputpath = outputdir + "sampled_graph.txt"
    if args.algorithm == "answer":
        inputpath = args.datasetdir + args.dataset + "_final.txt"
    assert(os.path.isfile(inputpath)), inputpath

    runflag = True
    if os.path.isfile(outputdir + "pairdeg_avg.txt"):
        if args.recalculate:
            os.remove(outputdir + "pairdeg_avg.txt")
        else:
            runflag = False
            num_line = 0
            with open(outputdir + "pairdeg_avg.txt", "r") as f:
                for line in f.readlines():
                    h, _ = line.rstrip().split(",")
                    if num_line >= len(hyperedge_steps):
                        runflag = True
                        os.remove(outputdir + "pairdeg_avg.txt")
                        break
                    if int(h) != hyperedge_steps[num_line]:
                        runflag = True
                        os.remove(outputdir + "pairdeg_avg.txt")
                        break
                    num_line += 1
            if num_line < len(hyperedge_steps):
                runflag = True
                os.remove(outputdir + "pairdeg_avg.txt")
            elif num_line == len(hyperedge_steps):
                print("Exist ", outputdir + "pairdeg_avg.txt")
    if runflag:
        tmpdict = defaultdict(int)
        number_of_edges = 0
        step = 0
        
        start_time = time.time()
        if args.verbose:
            print("Calculate Pairdeg Avg")

        with open(inputpath, "r") as f:
            for line in f.readlines():
                line = line[:-1]
                hyperedge = [int(_n) for _n in line.split(",")]
                number_of_edges += 1
                hsize = len(hyperedge)
                for i in range(hsize):
                    for j in range(i+1, hsize):
                        ni, nj = hyperedge[i], hyperedge[j]
                        if ni < nj:
                            key = str(ni) + "_" + str(nj)
                        else:
                            key = str(nj) + "_" + str(ni)
                        tmpdict[key] += 1
                
                if number_of_edges == hyperedge_steps[step]:
                    if args.verbose and step % 10 == 0:
                        runtime = (time.time() - start_time) // 60
                        print("[" + str(step) + "] Time(min.) = ", runtime)
                    step += 1
                    total = len(tmpdict.keys()) # number of connected node pairs
                    avgpairdegree = 0
                    for k, v in tmpdict.items():
                        avgpairdegree += (v / total) 
                    with open(outputdir + "pairdeg_avg.txt", "+a") as f:
                        f.write(str(number_of_edges) + "," + str(avgpairdegree) + "\n")
        print("End Pairdegree Avg")
        