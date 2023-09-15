#include <random>
#include <iterator>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <cstdlib> // for std::and() and std::srand()
# include "hset.hpp"

#ifndef ALGORITHMTIHS_HPP
#define ALGORITHMTIHS_HPP
using namespace std;

class Algorithm_TIHS{
public:
    HyperGraph *graph;
    string outputdir;
    vector<int> pool;

    vector<bool> hedge_check;
    vector<bool> check_node;
    vector<int> check_hyperedge;

    Algorithm_TIHS(string outputdir, HyperGraph *graph){
        this->outputdir = outputdir;
        this->graph = graph;
        this->hedge_check.resize(graph->number_of_hedges);
        this->check_node.resize(graph->number_of_nodes);
        this->check_hyperedge.resize(graph->number_of_hedges);
    }
    HSet* run(double accuracy);
    int sample_hedge(vector<bool> hedge_check, int remain);
};
#endif