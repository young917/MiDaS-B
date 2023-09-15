#include <cmath>
#include <cassert>
#include <chrono>
#include <random>
#include <iterator>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <cstdlib> 
#include <queue>
#include "graph.hpp"

#ifndef HELPER_HPP
#define HELPER_HPP

using namespace std;

// clusteringcoef, densification, intersection, sizewcc, density, overlapness, degree_avg
class Helper {
public:
    HyperGraph *graph;
    vector<bool> hypergraph_masking;
    vector<bool> node_masking;
    vector<int> check;
    vector<bool> check_node;
    vector<int> ancestor;
    vector<int> group_size;
    int number_of_hedges;
    int number_of_nodes;
    int sum_of_hsizes;

    bool ist_flag = false;
    bool ist_size_flag = false;
    bool dsf_flag = false;
    bool cc_flag = false;
    bool wcc_flag = false;
    bool cc_approx = false;
    string dataname;

    long long intersect = 0;
    unordered_map<int,long long> intersection_size;
    // unordered_map<int,long long> degree; // degree -> frequent
    // unordered_map<string, bool> neighbor;
    set<string> neighbor;
    string outputdir;
    queue<int> visited;
    vector<int> nodes;

    vector<int> vec;
    vector<int>::iterator it;
    vector<string> tmp;

    Helper(set<int> subhypergraph, HyperGraph *graph, string outputdir, string algo_opt);
    ~Helper(){
        hypergraph_masking.clear();
        node_masking.clear();
        check.clear();
        check_node.clear();
        // degree.clear();
        neighbor.clear();
    }
    void get_intersection(void);
    void get_intersection_size(void);
    void get_densification(void);
    void get_average_clustering_coef(void);
    void get_average_clustering_coef_appx(void);
    // void get_degreedist_property(void);
    void count_wcc(void);
    void get_dense_property(void);

    void update(set<int> deltaset, HyperGraph *graph);
    void save_properties(void);

private:
    int bfs(int start_node);
    int get_ancestor(int node);
};
#endif