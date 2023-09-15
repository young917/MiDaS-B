#include <random>
#include <iterator>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <cstdlib> // for std::and() and std::srand()

# include "hset.hpp"

#ifndef ALGORITHMES_HPP
#define ALGORITHMES_HPP

using namespace std;

class AlgorithmES {
public:
    int turn;
    double alpha;
    string outputdir;
    string algo_opt;
    int algo_opt_length;
    bool addflag;
    HyperGraph *graph;
    double resolution = 1.0;

    vector<int> pool;

    vector<int> hedge2prop; //current = (min / max/ avg), homogeneity
    // vector<int> node2prop; // current = (global deg)
    vector<int> hedges;
    vector<double> hedge_tree;
    int leaf_start;
    int max_key;
    int index = 0;
    unordered_map<int, vector<int>> leaf2hedges;
    // unordered_map<string, int> pairdegree;

    AlgorithmES(string outputdir, string algo_opt, double alpha, HyperGraph *graph){
        this->outputdir = outputdir;
        vector<string> options = split(algo_opt, '_');
        if (options[0].compare(0, 3, "add") == 0){
            this->addflag = true; 
            this->alpha = alpha;
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
        this->algo_opt_length = (int)this->algo_opt.size();
        this->graph = graph;
        
        // this->node2prop.resize(graph->number_of_nodes);
        this->hedge2prop.resize(graph->number_of_hedges);
        
        string parampath = outputdir + "parameters.txt";
        ofstream paramFile(parampath.c_str());
        paramFile << "algo opt: " << algo_opt << endl;
        paramFile << "alpha: " << to_string(alpha) << endl;
        paramFile.close();
    }
    ~AlgorithmES(){
        for (auto &item: leaf2hedges){
            auto &vec = item.second;
            vec.clear();
        }
        leaf2hedges.clear();
        hedge_tree.clear();
        hedge2prop.clear();
        // node2prop.clear();
        // hedges.clear();
    }

    HSet* run(double accuracy);
    void initiate(void);
    void update_prop(int &hprop, int cand);
    int sample_hedge(void);
};
#endif