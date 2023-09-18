# MiDaS-B: Back-In-Time Sampling from Real-world Hypergraphs


We extends our previous paper: [MiDaS: Representative Sampling from Real-world Hypergraphs](https://arxiv.org/abs/2202.01587) for back-in-time hypergraph sampling. 
*Back-in-time hypergraph sampling* is to accurately approximate a past snapshot of a given size of the input large hypergraph. 
Unlike representative sampling, we do not have access to the target (i.e., the past snapshot of a given size). 
To tackle this challenge, we propose **MiDaS-B**, an extension of MiDaS specifically designed for back-in-time hypergraph sampling. 

We provide source code for (1) sampling hypergraphs, (2) evaluating the quality of sub-hypergraphs for this back-in-time hypergraph sampling problem, and (3) visualizing results.


(1) **Sampling Hypergraphs**

* *Simple and Intuitive Approaches* : seven intuitive approaches(RNS, RDN, RW, FF, ONS, RHS, and TIHS) experimenting on 11 datasets
* *MiDaS* : a sampling method initially proposed for representative hypergraph sampling which allows for control over biases towards high-degree nodes through the parameter alpha
* *MiDaS-B*: an extension of MiDaS designed explicitly for back-in-time hypergraph sampling which incorporates a hyperedge-size-related term into the hyperedge sampling probabilities, effectively controlling associated biases

(2) **Evaluating the Quality of Sub-Hypergraph**

Measure how precisely the sub-hypergraph preserves the structural properties of the ground-truth past snapshot of the original hypergraph with respect to,

* *node-level statistics*: the distributions of node degrees and node-pair degrees
* *hyperedge-level statistics* : the distributions of hyperedge sizes and intersection sizes
* *graph-level statistics* : average clustering coefficient, density, overlapness, and effective diameter


(3) **Result Visualization**

Provide how to generate figures that illustrate 10 properties of the sub-hypergraphs obtained from different sampling approaches

- - -

## Datasets

In the paper, we used datasets after removing duplicated hyperedges. We preprocessed eleven datasets collected by [Austin R. Benson](https://www.cs.cornell.edu/~arb/data/). The datasets used in the paper are available in the ["dataset" folder](https://github.com/young917/MiDaS/tree/main/dataset)

- - -

## How to Run

Before you proceed with running the code, please ensure that you have compiled it using the `make` command

You can execute the following scripts:

* **All Baseline Sampling Methods**: Use [run_baseline.sh](https://github.com/young917/MiDaS-B/blob/main/run/run_baseline.sh) to run all seven baseline sampling methods

* **MiDaS w/o. Oracle**: Execute [run_midas.sh](https://github.com/young917/MiDaS-B/blob/main/run/run_midas.sh) to run MiDaS without utilizing an Oracle

* **MiDaS w. Oracle**: Utilize [run_midas_oracle.sh](https://github.com/young917/MiDaS-B/blob/main/run/run_midas_oracle.sh) to run MiDaS with Oracle support

* **MiDaS-B**: For our proposed method, MiDaS-B, execute [run_midasB.sh](https://github.com/young917/MiDaS-B/blob/main/run/run_midasB.sh)

* **MiDaS-B w. Oracle**: If you wish to run MiDaS-B with Oracle, use [run_midasB_oracle.sh](https://github.com/young917/MiDaS-B/blob/main/run/run_midasB_oracle.sh)

These scripts offer a convenient way to execute the sampling code and evaluate the output at once
 
### Evaluation

Before evaluating sampling methods, compute 10 properties of the ground-truth past snapshot of the original hypergraph by following [run_answer_analysis.sh](https://github.com/young917/MiDaS-B/blob/main/run/run_answer_analysis.sh)

```
# Density, CC, GCC
./bin/Sampling --dataname [dataname] --inputpath ../dataset/ --algorithm helper --inputdir [target_algo] --algo_opt "intersection,densification,sizewcc,clusteringcoef" --accuracy 100 --repeat [1,2,3]

# Degree, Int.Size, Pair Degree, Size Dist.
./bin/Sampling --dataname [dataname] --inputpath ../dataset/ --algorithm helperdist --inputdir [target_algo]  --samplingportion [0.1,0.3,0.5,0.7,0.9] --repeat [1,2,3]

# Diameter, Overlapness
cd src
python calculation_helper.py --dataset [dataname] --algorithm [target_algo] --effdiam --overlapness --repeat [1,2,3] --accuracy 100

# SV
# For threads-ask-ubuntu, coauth-MAG-Geology-full, and coauth-MAG-History-full, 
# run  `src/preprocess_sv.py --dataname [dataname] --algoname [target_algo] --repeat_str [1,2,3]` 
# and then run `run_sv_midas_oracle.m` before the next line

cd src
python calculation_helper.py --dataset [dataname] --algorithm [target_algo] --svdist --samplingportion [0.1,0.3,0.5,0.7,0.9] --repeat [1,2,3]

```

After computing all 10 properties of the sampling results, run `analyze/evaluation_table.py --search_name [baseline,final,ablation] --repeat "1,2,3"`


### Visualization

You can replicate the table displaying the performance of various sampling methods across 11 datasets, considering 10 properties and 5 distinct sampling portions. 
To do so, you can refer to the [Performance.ipynb](https://github.com/young917/MiDaS-B/blob/main/analyze/Performance.ipynb)

Furthermore, engage in a detailed comparison of the properties exhibited by sub-hypergraphs generated by different sampling methods. This comparison is available in the [PlotDistribution.ipynb](https://github.com/young917/MiDaS-B/blob/main/analyze/PlotDistribution.ipynb)

