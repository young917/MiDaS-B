#include <random>
#include <iterator>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <cstdlib> // for std::and() and std::srand()

# include "hset.hpp"

#ifndef ALGORITHMESSZ_HPP
#define ALGORITHMESSZ_HPP

using namespace std;

class AlgorithmESSZ {
public:
    double alpha;
    double beta;
    string outputdir;
    string algo_opt;
    int algo_opt_length;
    bool addflag;
    HyperGraph *graph;
    double resolution = 1.0;

    vector<int> pool;

    vector<int> hedge2prop; //current = (min / max/ avg), homogeneity
    unordered_map<string, int> prop2index;
    unordered_map<int, string> index2prop;
    vector<int> node2prop; // current = (global deg)
    // vector<double> prob_dist;
    vector<double> hedge_tree;
    int leaf_start;
    int max_key;
    unordered_map<int, vector<int>> leaf2hedges;

    AlgorithmESSZ(string outputdir, string algo_opt, double alpha, double beta, HyperGraph *graph){
        this->outputdir = outputdir;
        vector<string> options = split(algo_opt, '_');
        if (options[0].compare(0, 3, "add") == 0){
            this->addflag = true; 
            this->alpha = alpha;
            this->beta = -beta;
        }
        else if (options[0].compare(0, 3, "rem") == 0){
            this->addflag = false;
            this->alpha = -alpha;
            this->beta = beta;
        }
        int len = (int)options.size();
        this->algo_opt = options[1];
        for (int i = 2; i < len ; i++){
            this->algo_opt += "_";
            this->algo_opt += options[i];
        }
        this->algo_opt_length = (int)this->algo_opt.size();
        this->graph = graph;
        
        this->node2prop.resize(graph->number_of_nodes);
        this->hedge2prop.resize(graph->number_of_hedges);
        
        string parampath = outputdir + "parameters.txt";
        ofstream paramFile(parampath.c_str());
        paramFile << "algo opt: " << algo_opt << endl;
        paramFile << "alpha: " << to_string(alpha) << "  beta: " << to_string(beta) << endl;
        paramFile.close();
    }
    ~AlgorithmESSZ(){
        for (auto &item: leaf2hedges){
            auto &vec = item.second;
            vec.clear();
        }
        leaf2hedges.clear();
        hedge_tree.clear();
        hedge2prop.clear();
        node2prop.clear();
        prop2index.clear();
        index2prop.clear();
    }

    HSet* run(double accuracy);
    void initiate(void);
    void update_prop(int &hprop, int cand);
    int sample_hedge(void);

};
#endif