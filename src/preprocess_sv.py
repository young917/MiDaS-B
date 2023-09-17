from collections import defaultdict
from tqdm import tqdm
import math
import os
import shutil
import argparse

def find_incidence_mat(inputpath, outputname, portion):
    node2edges = defaultdict(list)
    nodename2index = {}
    node_index = 0
    number_of_nodes = 0
    number_of_edges = 0
    if os.path.isfile(inputpath) is False:
        with open("not_exist.txt", "+a") as f:
            f.write(inputpath + "\n")
        return -1
         
    with open(inputpath, "r") as f:
        for idx, line in enumerate(f.readlines()):
            line = line[:-1] # strip enter
            nodes = line.split(",")
            for i, node in enumerate(nodes):
                if node not in nodename2index:
                    nodename2index[node] = node_index
                    node_index += 1
                nodes[i] = nodename2index[node]
            for v in nodes:
                node2edges[int(v)].append(idx)
            number_of_edges += 1
    number_of_nodes = len(node2edges.keys())

    if number_of_edges == 0:
        print("Empty Input")
        return -1

    sampled_numhyperedge = int(math.ceil(number_of_edges * portion))
    print(number_of_nodes, number_of_edges)
    total = 0
    with open(outputname + "dim_%.1f.txt" % (portion), "w") as f:
        f.write(str(number_of_nodes) + "\n")
        f.write(str(sampled_numhyperedge) + "\n")
        # f.write(str(number_of_edges) + "\n")
    with open(outputname + "icd_row_%.1f.txt" % (portion), "w") as fr, open(outputname + "icd_col_%.1f.txt" % (portion), "w") as fc:
        for v in tqdm(range(number_of_nodes)):
            for h in node2edges[v]:
                if h < sampled_numhyperedge:
                    fr.write(str(v + 1) + "\n")
                    fc.write(str(h + 1) + "\n")
                    total += 1
    with open(outputname + "icd_%.1f.txt" % (portion), "w") as f:
        for v in tqdm(range(number_of_nodes)):
            for h in node2edges[v]:
                if h < sampled_numhyperedge:
                    f.write(str(v+1) + ":" + str(h+1) + "\n")
                    total += 1
    print(total)
            
if __name__ == "__main__":
    dataset = ["threads-ask-ubuntu", "coauth-MAG-Geology-full", "coauth-MAG-History-full"]

    sampled_algo_list = ["ns/add_global_deg_0.0000", "ns/add_global_deg_1.0000", "rw/rw_c_1.00_0.30", "ff/ff_c_0.51_0.20", 'es/add_global_deg_min_0.0000', "ns/add_ordered", "tihs"]
    
    # For midas - oracle
    # alphas = ["0.0000", "0.1250", "0.2500", "0.5000", "1.0000", "2.0000"]
    # sampled_alpha_algo_list = ["es/add_global_deg_min_"] 
    # for alpha in alphas:
    #     for algo in sampled_alpha_algo_list:
    #         sampled_algo_list.append(algo + alpha)

    # For midasB - oracle
    # alphas = ["0.00", "0.12", "0.25", "0.50", "1.00", "2.00"]
    # betas = ["-1.00", "-0.50", "-0.25", "-0.12", "0.00", "0.12", "0.25", "0.50", "1.00"]
    # sampled_alpha_beta_algo_list = ["essz/add_global_deg_min_"]
    # for alpha in alphas:
    #     for beta in betas:
    #         for algo in sampled_alpha_beta_algo_list:
    #             sampled_algo_list.append(algo + alpha + "_" + beta)
    

    parser = argparse.ArgumentParser()
    parser.add_argument('--repeat_str', required=False, default="1,2,3")
    parser.add_argument('--portion_str', required=False, default="0.1,0.3,0.50,0.70,0.90")
    parser.add_argument('--dataname', required=False, default="")
    parser.add_argument('--algoname', required=False, default="")
    parser.add_argument('--recalculate', action='store_true')

    args = parser.parse_args()
    repeat_list = args.repeat_str.split(",")
    portion_list = args.portion_str.split(",")
    if len(args.dataname) > 0:
        dataset = args.dataname.split(",")
    if len(args.algoname) > 0:
        sampled_algo_list = args.algoname.split(",")
    
    for dataname in dataset:
        for portion_str in portion_list:
            for algoname in tqdm(sampled_algo_list, desc=dataname + " " + portion_str):
                portion = float(portion_str)
                if algoname == "answer":
                    inputpath = "../../dataset/" + dataname + ".txt"
                    outputdir = "../analyze_sv/input/answer/" + dataname + "/"
                    
                    if os.path.isdir(outputdir) is False:
                        os.makedirs(outputdir)
                    
                    if (args.recalculate is False) and (os.path.isfile("../results/answer/" + dataname + "/singular_values_full_%.1f.txt" % (portion))):
                        #(os.path.isfile("../Hypergraph_Sampling_cpp/results/answer_dist/" + dataname + "/singular_values_full.txt")):
                        shutil.rmtree(outputdir)
                    else:    
                        later_outputdir = "../analyze_sv/output/answer/" + dataname + "/"
                        if os.path.isdir(later_outputdir) is False:
                            os.makedirs(later_outputdir)
                        answer_s = find_incidence_mat(inputpath, outputdir, portion)
                else:                    
                    for repeat_str in repeat_list:
                        inputpath = "../results/" + algoname + "/" + dataname + "/" + repeat_str + "/sampled_graph.txt"
                        outputdir = "../analyze_sv/input/" + algoname + "/" + dataname  + "/" + repeat_str + "/"
                        if os.path.isdir(outputdir) is False:
                            os.makedirs(outputdir)
                        if (args.recalculate is False) and (os.path.isfile("../analyze_sv/output/" + algoname + "/" + dataname + "/" + repeat_str + "/singular_values_full_%.1f.txt" % (portion))):
                            #(os.path.isfile("../Hypergraph_Sampling_cpp/results/" + algoname + "/" + dataname + "_" + portion_str  + "/" + repeat_str + "/singular_values_full.txt")):
                            shutil.rmtree(outputdir)
                        else:    
                            later_outputdir = "../analyze_sv/output/" + algoname + "/" + dataname + "/" + repeat_str + "/"
                            if os.path.isdir(later_outputdir) is False:
                                os.makedirs(later_outputdir)
                            data_s = find_incidence_mat(inputpath, outputdir, portion)
                            