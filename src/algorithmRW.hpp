# include <random>
# include <queue>
# include "hset.hpp"

#ifndef ALGORITHMRW_HPP
#define ALGORITHMRW_HPP
using namespace std;

class Algorithm_RW{
public:
    HyperGraph *graph;
    double restart;
    string algo_opt, outputdir;
    double givenmaxlength;
    bool noinduce;
    vector<vector<int>> node2node;
    vector<int> htable;
    vector<bool> check;

    vector<int> pool;
    queue<int> q;
    vector<int> order;
    vector<int> tmp;

    Algorithm_RW(string algo_opt, string outputdir, HyperGraph *graph, double restart, double givenmaxlength, bool noinduce);
    HSet* run(double accruacy);
    void walk(int seed_node, int max_length);
};
#endif