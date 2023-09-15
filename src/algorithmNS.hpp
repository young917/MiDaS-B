#include <random>
#include <iterator>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <cstdlib> // for std::and() and std::srand()

# include "hset.hpp"

#ifndef ALGORITHMNS_HPP
#define ALGORITHMNS_HPP

using namespace std;

class AlgorithmNS {
public:
    string outputdir;
    string algo_opt;
    HyperGraph *graph;
    double alpha;
    bool addflag;

    vector<int> nodes;
    vector<double> node_tree;
    int leaf_start;
    unordered_map<int, vector<int>> leaf2nodes;

    vector<int> htable;
    vector<int> pool;

    int index = 0;
    int order = 0;

    AlgorithmNS(double alpha, string outputdir, string algo_opt, HyperGraph *graph){
        this->alpha = alpha;
        this->outputdir = outputdir;
        this->graph = graph;
        this->htable.resize(graph->number_of_hedges); // For induced hyperedges)
        vector<string> options = split(algo_opt, '_');
        if (options[0].compare(0, 3, "add") == 0){
            this->addflag = true; 
        }
        else if (options[0].compare(0, 3, "rem") == 0){
            this->addflag = false;
            this->alpha = -alpha;
        }
        int len = (int)options.size();
        this->algo_opt = options[1];
        for (int i = 2; i < len ; i++){
            this->algo_opt += "_";
            this->algo_opt += options[i];
        }
    }
    HSet* run(double accuracy);

    int sample_he(HSet *sampled, int flag);
    void initiate(void);
    void reinitiate(void);
    int sample_node(void);
};
#endif