import os
import sys
import argparse
from collections import defaultdict
from tqdm import tqdm
import numpy as np
from itertools import chain
import networkx as nx
from numpy.linalg import svd
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import svds, eigs, norm
from scipy.stats import skew, kendalltau
import snap
import powerlaw
import time
import math
import copy

# ----------------------------------------------------------------------------------------------
num = 300

def find_svs(dataname, inputpath, outputdir, recalculate_flag, samplingportion, numhyperedge):
    outputpath = outputdir + "sv_full_%.1f.txt" % (samplingportion)
    outputpath2 = '../analyze_sv/output/' + outputdir[11:] + "sv_full_%.1f.txt" % (samplingportion)
    node2edges = defaultdict(list)
    nodename2index = {}
    node_index = 0
    number_of_nodes = 0
    number_of_edges = 0

    print(inputpath)

    with open(inputpath, "r") as f:
        for idx, line in enumerate(f.readlines()):
            line = line[:-1] # strip enter
            nodes = line.split(",")
            for i in range(len(nodes)):
                node = nodes[i]
                if node not in nodename2index:
                    nodename2index[node] = node_index
                    node_index += 1
                nodes[i] = nodename2index[node]
            for v in nodes:
                node2edges[int(v)].append(idx)
            number_of_edges += 1
            if number_of_edges == numhyperedge:
                break
    number_of_nodes = len(node2edges.keys())

    if number_of_edges == 0:
        print("Empty Input")
        return -1, -1, -1

    dim = min(number_of_edges, number_of_nodes)
    s = -1
    flag = True
    if (recalculate_flag is False):
        if (os.path.isfile(outputpath)): # Read Singular Values
            s = []
            with open(outputpath, "r") as f:
                for line in f.readlines():
                    line = line[:-1]
                    s.append(float(line))
                    
            # if (dataname in exceptdatas) is False and len(s) < dim: # Error -> Recalculate
            #     flag = True
            if len(s) < dim:
                flag = True
            else:
                flag = False
        elif (os.path.isfile(outputpath2)):
            s = []
            with open(outputpath2, "r") as f:
                for line in f.readlines():
                    line = line[:-1]
                    s.append(float(line))
            assert len(s) == 300, len(s)
            flag = False
            dim = 300
    # elif (recalculate_flag is False) and (dataname in exceptdatas):
    #     print("no singular value file", outputpath)
    #     flag = False
        
    # Calculate Singular Values
    if (flag): 
        try:
            incident_matrix = np.zeros(shape=(number_of_nodes, number_of_edges), dtype=np.byte)
            for v in node2edges.keys():
                for edge_idx in node2edges[v]:
                    incident_matrix[v,edge_idx] = 1
            # sum_of_squares = LA.norm(incident_matrix, 'fro') ** 2            
            s = svd(incident_matrix, compute_uv=False)
            s = sorted(list(s), reverse=True)
            assert len(s) == dim
            for i in range(1,dim):
                assert s[i-1] >= s[i]
            # assert math.fabs(sum_of_squares - np.sum(np.array(s) ** 2)) < 0.000001
            with open(outputpath, "w") as f:
                for _sv in s:
                    f.write(str(_sv) + "\n")
        except:
            # print("[" + inputpath + "] Error #V=" + str(number_of_nodes) + ", #E=" + str(number_of_edges))
            rows, cols = zip(*chain.from_iterable([[(v, edge_idx) for edge_idx in node2edges[v]] for v in node2edges.keys()]))
            nnz = len(rows)
            incident_matrix = coo_matrix((np.ones(nnz), (rows, cols)), shape=(number_of_nodes, number_of_edges))
            sum_of_squares = norm(incident_matrix, 'fro') ** 2
            rank = min(num , dim - 1)
            _, s, _ = svds(incident_matrix.tocsc(), k=rank)
            last_sv_square = sum_of_squares - sum([_s * _s for _s in s])
            last_sv = math.sqrt(last_sv_square)
            s = list(s)
            s.append(last_sv)
            assert len(s) == dim
            s = sorted(s, reverse=True)
            for i in range(1,dim):
                assert s[i-1] >= s[i]
            with open(outputpath, "w") as f:
                for _sv in s:
                    f.write(str(_sv) + "\n")

    return s, dim

def sv_dist(_list_sv, dim, answer_max_portion):
    list_sv = copy.deepcopy(_list_sv)
    singular_values = {}
    
    if answer_max_portion != -1:
        number_of_required_svs = math.ceil(dim * answer_max_portion)
        list_sv = list_sv[:number_of_required_svs]
    denom = sum([_s * _s for _s in list_sv])
    until = 0.0
    total = dim
    for idx in range(len(list_sv)):
        proportion = (idx + 1) / total
        until += list_sv[idx] ** 2
        singular_values[proportion] = until / denom

    return singular_values

# ----------------------------------------------------------------------------------------------

RESULTDIR = "/home/MiDaS-B/results/"

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', required=False, default="email-Enron-full")
parser.add_argument('--datasetdir', required=False, default="../../dataset/")
parser.add_argument('--algorithm', required=False, default='es/global_deg_min_0.0000')
parser.add_argument('--repeat', required=False, type=int, default=-1)
parser.add_argument('--accuracy', required=False, type=int, default=500)


parser.add_argument('--effdiam', required=False, action='store_true')
parser.add_argument('--sv', required=False, action='store_true')
parser.add_argument('--svdist', required=False, action='store_true')
parser.add_argument('--overlapness', required=False, action='store_true')

parser.add_argument('--verbose', required=False, action='store_true')
parser.add_argument('--recalculate', required=False, action='store_true')

parser.add_argument('--samplingportion', required=False, type=float, default=0.5)

args = parser.parse_args()
outputdir = '../results/' + args.algorithm + "/" + args.dataset + "/"
if args.repeat > 0:
    outputdir += str(args.repeat) + "/"
print("Start " + outputdir)

appx_flag = False
appx_flag2 = False
if args.dataset in ["threads-ask-ubuntu", "tags-ask-ubuntu", "tags-math-sx"]:
    appx_flag = True
elif args.dataset in ["coauth-DBLP-full", "threads-math-sx", "coauth-MAG-Geology-full"]:
    appx_flag2 = True

entire_nodeset = set([])
entire_number_of_hyperedges = 0
datapath = args.datasetdir + args.dataset + "_final.txt"
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

if args.effdiam:
    hyperedge_steps = []
    path = RESULTDIR + "answer/" + args.dataset + "/effdiameter.txt"
    if os.path.isfile(path) is False:
        path = RESULTDIR + "answer/" + args.dataset + "/clusteringcoef.txt"
        with open(path, "r") as f:
            for line in f.readlines():
                numh = int(line.rstrip().split(",")[0])
                hyperedge_steps.append(numh)
    else:
        with open(path, "r") as f:
            for line in f.readlines():
                numh = int(line.rstrip().split(",")[1])
                hyperedge_steps.append(numh)

    inputpath = outputdir + "sampled_graph.txt"
    if args.algorithm == "answer":
        inputpath = datapath
    step = 0 #thr = args.accuracy
    assert(os.path.isfile(inputpath)), inputpath

    runflag = True
    if os.path.isfile(outputdir + "effdiameter.txt"):
        if args.recalculate:
            os.remove(outputdir + "effdiameter.txt")
        else:
            runflag = False
            num_line = 0
            with open(outputdir + "effdiameter.txt", "r") as f:
                for line in f.readlines():
                    _, h, _ = line.rstrip().split(",")
                    if num_line >= len(hyperedge_steps):
                        runflag = True
                        os.remove(outputdir + "effdiameter.txt")
                        break
                    if int(h) != hyperedge_steps[num_line]:
                        runflag = True
                        os.remove(outputdir + "effdiameter.txt")
                        break
                    num_line += 1
            if num_line < len(hyperedge_steps) and runflag is False:
                runflag = True
                os.remove(outputdir + "effdiameter.txt")
            elif num_line == len(hyperedge_steps):
                print("Exist ", outputdir + "effdiameter.txt")

    if runflag:
        step = 0
        # build snap graph
        pg = snap.TUNGraph.New()
        nodeset = set([])
        number_of_edges = 0

        start_time = time.time()
        if args.verbose:
            print("Calculate Effective Diameter")
        with open(inputpath, "r") as f:
            for line in f.readlines():
                line = line[:-1]
                hyperedge = [int(_n) for _n in line.split(",")]
                number_of_edges += 1
                # Reindex
                # for i, node in enumerate(hyperedge):
                #     if node not in nodename2index:
                #         nodename2index[node] = node_index
                #         node_index += 1
                #     hyperedge[i] = nodename2index[node]
                # Add node
                for i, n in enumerate(hyperedge):
                    nodeindex = n
                    if args.algorithm == "answer":
                        nodeindex = datanode2index[n]
                        hyperedge[i] = nodeindex
                    if nodeindex not in nodeset:
                        nodeset.add(nodeindex)
                        pg.AddNode(nodeindex)
                # Add clique edges
                if len(hyperedge) != 1: 
                    for i in range(0, len(hyperedge)-1):
                        for j in range(i+1, len(hyperedge)):
                            i1, i2 = min(hyperedge[i], hyperedge[j]), max(hyperedge[i], hyperedge[j])
                            ret = pg.AddEdge(i1, i2)
                num_nodes = len(nodeset)
                # if (num_nodes > 10) and ((num_nodes / entire_number_of_nodes) > thr):
                #if ((number_of_edges / entire_number_of_hyperedges) >= thr) or (number_of_edges == entire_number_of_hyperedges):
                if (number_of_edges == hyperedge_steps[step]):
                    if args.verbose and step % 10 == 0:
                        runtime = (time.time() - start_time) // 60
                        print("[" + str(step) + "], ", num_nodes, " Time(min) : ", runtime)
                    step += 1 # thr += args.accuracy
                    if appx_flag:
                        effective_diameter = snap.GetBfsEffDiam(pg, num_nodes if num_nodes < 5000 else 1000, False)
                    elif appx_flag2:
                        if step <= 5:
                            effective_diameter = snap.GetBfsEffDiam(pg, num_nodes if num_nodes < 5000 else 1000, False)
                        else:
                            effective_diameter = snap.GetBfsEffDiam(pg, num_nodes if num_nodes < 5000 else 100, False)
                    else:
                        effective_diameter = snap.GetBfsEffDiam(pg, num_nodes if num_nodes < 5000 else 5000, False)
                    # Save Effective Diameter
                    with open(outputdir + "effdiameter.txt", "a+") as f:
                        f.write(str(num_nodes) + "," + str(number_of_edges) + "," + str(effective_diameter) + "\n")
        print("End EffDiameter")

if args.sv:
    hyperedge_steps = []
    unit_numhedge = int(math.ceil(entire_number_of_hyperedges / args.accuracy))
    for _i in range(entire_number_of_hyperedges):
        if (_i > 0) and (_i % unit_numhedge) == 0:
            hyperedge_steps.append(_i)
    hyperedge_steps.append(entire_number_of_hyperedges)
                
    inputpath = outputdir + "sampled_graph.txt"
    if args.algorithm == "answer":
        inputpath = datapath
    assert(os.path.isfile(inputpath))
    runflag = True
    if os.path.isfile(outputdir + "sv.txt"):
        if args.recalculate:
            os.remove(outputdir + "sv.txt")
        else:
            runflag = False
            num_line = 0
            with open(outputdir + "sv.txt", "r") as f:
                for line in f.readlines():
                    _, h, _ = line.rstrip().split(",")
                    if num_line >= len(hyperedge_steps):
                        runflag = True
                        os.remove(outputdir + "sv.txt")
                        break
                    if int(h) != hyperedge_steps[num_line]:
                        runflag = True
                        os.remove(outputdir + "sv.txt")
                        break
                    num_line += 1
            if num_line < len(hyperedge_steps) and runflag is False:
                runflag = True
                os.remove(outputdir + "sv.txt")
            elif num_line == len(hyperedge_steps):
                print("Exist ", outputdir + "sv.txt")
    if runflag:
        step = 0 #thr = args.accuracy
        node2edge = defaultdict(list)
        oldindex2new = {}
        node_index = 0
        number_of_nodes = 0
        number_of_edges = 0

        start_time = time.time()
        if args.verbose:
            print("Calculate SV")
        with open(inputpath, "r") as f:
            for hidx, line in enumerate(f.readlines()):
                line = line[:-1]
                nodes = [int(_n) for _n in line.split(",")]
                # Reindex
                for i, node in enumerate(nodes):
                    if node not in oldindex2new:
                        oldindex2new[node] = node_index
                        node_index += 1
                    nodes[i] = oldindex2new[node]
                for node in nodes:
                    node2edge[node].append(hidx)
                number_of_edges += 1
                number_of_nodes = node_index
                
                # if ((number_of_edges / entire_number_of_hyperedges) >= thr) or (number_of_edges == entire_number_of_hyperedges):
                if hyperedge_steps[step] == number_of_edges:
                    if args.verbose and step % 10 == 0:
                        runtime = (time.time() - start_time) // 60
                        print("[" + str(step) + "], Time(min) : ", runtime)
                    step += 1 # thr += args.accuracy
                    rows, cols = zip(*chain.from_iterable([[(v, edge_idx) for edge_idx in edges]
                                                            for v, edges in node2edge.items()]))
                    nnz = len(rows)
                    incident_matrix = coo_matrix((np.ones(nnz), (rows, cols)), shape=(number_of_nodes, number_of_edges))
                    # print(incident_matrix.shape)
                    _, s, _ = svds(incident_matrix.tocsc(), k=1)
                    # Save Largest Singular Values
                    with open(outputdir + "sv.txt", "a+") as f:
                        f.write(str(number_of_nodes) + "," + str(number_of_edges) + "," + str(s[0]) + "\n")

        print("End Singular Values")

if args.svdist:
    numhyperedge = int(math.ceil(entire_number_of_hyperedges * args.samplingportion))
    inputpath = outputdir + "sampled_graph.txt"
    if args.algorithm == "answer":
        inputpath = datapath
        assert(os.path.isfile(inputpath))
        answer_s, answer_dim = find_svs(args.dataset, inputpath, outputdir, args.recalculate, args.samplingportion, numhyperedge)
        answer_dict = sv_dist(answer_s, answer_dim, -1)
        with open(outputdir + "svdist_%.1f.txt" % (args.samplingportion), "w") as f:
            for k,v in answer_dict.items():
                f.write(str(k) + " " + str(v) + "\n")
    else:
        assert(os.path.isfile(inputpath))
        sampled_s, sampled_dim = find_svs(args.dataset, inputpath, outputdir, args.recalculate, args.samplingportion, numhyperedge)             
        sampled_dict = sv_dist(sampled_s, sampled_dim, -1)
        with open(outputdir + "svdist_%.1f.txt" % (args.samplingportion), "w") as f:
            for k,v in sampled_dict.items():
                f.write(str(k) + " " + str(v) + "\n")

    print("End Singular Value Dist")

if args.overlapness:
    hyperedge_steps = []
    path = RESULTDIR + "answer/" + args.dataset + "/overlapness.txt"
    if os.path.isfile(path) is False:
        path = RESULTDIR + "answer/" + args.dataset + "/densification.txt"
        with open(path, "r") as f:
            for line in f.readlines():
                numh = int(line.rstrip().split(",")[-1])
                hyperedge_steps.append(numh)
    else:
        with open(path, "r") as f:
            for line in f.readlines():
                numh = int(line.rstrip().split(",")[0])
                hyperedge_steps.append(numh)

    inputpath = outputdir + "sampled_graph.txt"
    if args.algorithm == "answer":
        inputpath = datapath
    assert(os.path.isfile(inputpath))
    runflag = True
    if os.path.isfile(outputdir + "overlapness.txt"):
        if args.recalculate:
            os.remove(outputdir + "overlapness.txt")
        else:
            runflag = False
            num_line = 0
            with open(outputdir + "overlapness.txt", "r") as f:
                for line in f.readlines():
                    h, _ = line.rstrip().split(",")
                    if num_line >= len(hyperedge_steps):
                        runflag = True
                        os.remove(outputdir + "overlapness.txt")
                        break
                    if int(h) != hyperedge_steps[num_line]:
                        runflag = True
                        os.remove(outputdir + "overlapness.txt")
                        break
                    num_line += 1
            if num_line < len(hyperedge_steps) and runflag is False:
                runflag = True
                os.remove(outputdir + "overlapness.txt")
            elif num_line == len(hyperedge_steps):
                print("Exist ", outputdir + "overlapness.txt")
    if runflag:
        step = 0 # thr = args.accuracy
        oldindex2new = {}
        node_index = 0
        number_of_nodes = 0
        number_of_edges = 0
        sum_of_hyperedge_sizes = 0
        with open(inputpath, "r") as f:
            for hidx, line in enumerate(f.readlines()):
                line = line[:-1]
                nodes = [int(_n) for _n in line.split(",")]
                # Reindex
                for i, node in enumerate(nodes):
                    if node not in oldindex2new:
                        oldindex2new[node] = node_index
                        node_index += 1
                number_of_edges += 1
                number_of_nodes = node_index
                sum_of_hyperedge_sizes += len(nodes)
                
                # if ((number_of_edges / entire_number_of_hyperedges) >= thr) or (number_of_edges == entire_number_of_hyperedges):
                if hyperedge_steps[step] == number_of_edges:
                    if args.verbose and step % 10 == 0:
                        print("[" + str(step) + "]")
                    step += 1 # thr += args.accuracy
                    overlapness = sum_of_hyperedge_sizes / number_of_nodes
                    with open(outputdir + "overlapness.txt", "a+") as f:
                        f.write(str(number_of_edges) + "," + str(overlapness) + "\n")

        print("End Overlapness")
