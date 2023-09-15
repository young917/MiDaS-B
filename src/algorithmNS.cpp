#include "algorithmNS.hpp"

void AlgorithmNS::initiate(void){
    if (algo_opt.compare(0, 6, "random") == 0){
        nodes.clear();
        for(int v = 0 ; v < graph->number_of_nodes ; v++){
            nodes.push_back(v);
        }
    }
    else if (algo_opt.compare(0, 7, "ordered") == 0){
        nodes.clear();
        for(int v = 0 ; v < graph->number_of_nodes ; v++){
            nodes.push_back(v);
        }
    }
    else if (algo_opt.compare(0, 6, "decdeg") == 0){
        vector<pair<int,int>> tmp;
        for(int v = 0 ; v < graph->number_of_nodes ; v++){
            int vdeg = (int)graph->node2hyperedge[v].size();
            tmp.push_back(make_pair(-vdeg, v));
        }
        sort(tmp.begin(), tmp.end());
        nodes.clear();
        for(int vi = 0 ; vi < graph->number_of_nodes ; vi++){
            int v = tmp[vi].second;
            nodes.push_back(v);
        }
    }
    // dynamic node degree impossible..!
    else if(algo_opt.compare(0, 10, "global_deg") == 0){
        // clean data
        for(auto &item : leaf2nodes){
            auto &vec = item.second;
            vec.clear();
        }
        leaf2nodes.clear();
        node_tree.clear();
        
        // resize map
        int max_deg = 0;
        for(int v = 0; v < graph->number_of_nodes ; v++){
            int nodedegree = (int)graph->node2hyperedge[v].size();
            max_deg = max(max_deg, nodedegree);
        }
        int height = (int)ceil(log2(max_deg + 1));
        int whole_size = exp2(height+1);
        leaf_start = exp2(height);    
        node_tree.resize(whole_size, 0);
        
        // store
        for(int v = 0; v < graph->number_of_nodes ; v++){
            int nodedegree = graph->node2hyperedge[v].size();
            leaf2nodes[nodedegree].push_back(v);
        }
        for(int d = 0 ; d <= max_deg ; d++){
            if ((d == 0) || ((int)leaf2nodes[d].size() == 0)){
                node_tree[leaf_start + d] = 0.0;
            }
            else if (alpha == 0){
                node_tree[leaf_start + d] = (double)leaf2nodes[d].size();
            }
            else{
                node_tree[leaf_start + d] = pow((double)d, alpha) * (double)leaf2nodes[d].size();
            }
            assert (node_tree[leaf_start + d] >= 0.0);
        }
        for(int p = leaf_start - 1 ; p > 0 ; p--){
            node_tree[p] = node_tree[2 * p] + node_tree[2 * p + 1];
            assert (node_tree[p] >= 0.0);
        }
    }
}

void AlgorithmNS::reinitiate(void){
    // only global_deg
    if(algo_opt.compare(0, 10, "global_deg") == 0){
        // clean data
        for(auto &item : leaf2nodes){
            auto &vec = item.second;
            vec.clear();
        }
        // 
        int max_deg = 0;
        for(int v = 0; v < graph->number_of_nodes ; v++){
            int nodedegree = (int)graph->node2hyperedge[v].size();
            max_deg = max(max_deg, nodedegree);
        }
        // store
        for(int v = 0; v < graph->number_of_nodes ; v++){
            int nodedegree = graph->node2hyperedge[v].size();
            leaf2nodes[nodedegree].push_back(v);
        }
        for(int d = 0 ; d <= max_deg ; d++){
            if ((d == 0) || ((int)leaf2nodes[d].size() == 0)){
                node_tree[leaf_start + d] = 0.0;
            }
            else if (alpha == 0){
                node_tree[leaf_start + d] = (double)leaf2nodes[d].size();
            }
            else{
                node_tree[leaf_start + d] = pow((double)d, alpha) * (double)leaf2nodes[d].size();
            }
            assert (node_tree[leaf_start + d] >= 0.0);
        }
        for(int p = leaf_start - 1 ; p > 0 ; p--){
            node_tree[p] = node_tree[2 * p] + node_tree[2 * p + 1];
            assert (node_tree[p] >= 0.0);
        }
    }
}

int AlgorithmNS::sample_node(void){
    std::random_device rd;
    std::mt19937 gen(rd());

    if (algo_opt.compare(0, 6, "random") == 0){
        int numnodes = (int)nodes.size();
        std::uniform_int_distribution<> dist(0, numnodes-1);
        int random_index = dist(gen);
        int sampled_node = nodes[random_index];
        nodes.erase(nodes.begin() + random_index);
        return sampled_node;
    }
    else if (algo_opt.compare(0, 7, "ordered") == 0){
        return nodes[order++];
    }
    else if (algo_opt.compare(0, 6, "decdeg") == 0){
        int sampled_node = nodes[index];
        index += 1;
        return sampled_node;
    }
    else if (algo_opt.compare(0, 10, "global_deg") == 0){
        std::uniform_real_distribution<> dist(0, 1);
        int idx = 1;
        while (idx < leaf_start){
            assert(node_tree[idx] > 0);
            if (node_tree[2 * idx] == 0){
                idx = 2 * idx + 1;
                continue;
            }
            else if (node_tree[2 * idx + 1] == 0){
                idx = 2 * idx;
                continue;
            }
            double random_double = dist(gen);
            assert ((node_tree[2 * idx] + node_tree[2 * idx + 1]) > 0.0);
            double prob = (double)node_tree[2 * idx] / (node_tree[2 * idx] + node_tree[2 * idx + 1]);
            if (random_double < prob) idx = 2 * idx;
            else idx = 2 * idx + 1;
            assert (node_tree[idx] > 0.0);
        }
        int final_index = idx - leaf_start;

        int numnodes = (int)leaf2nodes[final_index].size();
        std::uniform_int_distribution<> dist2(0, numnodes-1);
        int random_index = dist2(gen);
        int sampled_node = leaf2nodes[final_index][random_index];

        // erase
        leaf2nodes[final_index].erase(leaf2nodes[final_index].begin() + random_index);
        if ((final_index == 0) || ((int)leaf2nodes[final_index].size() == 0)){
            node_tree[leaf_start + final_index] = 0.0;
        }
        else if (alpha == 0){
            node_tree[leaf_start + final_index] = (double)leaf2nodes[final_index].size();
        }
        else{
            node_tree[leaf_start + final_index] = pow((double)final_index, alpha) * (double)leaf2nodes[final_index].size();
        }
        assert (node_tree[leaf_start + final_index] >= 0.0);
        int p = leaf_start + final_index;
        while (p > 1){
            int parent = (int)(floor(p/2));
            node_tree[parent] = node_tree[2 * parent] + node_tree[2 * parent + 1];
            assert (node_tree[parent] >= 0.0);
            p = parent;
        }
        return sampled_node;
    }
    return -1;
}

HSet* AlgorithmNS::run(double accuracy){
    random_device rd;
    mt19937 g(rd());

    // Time Measure,
    // auto start = chrono::steady_clock::now();
    // auto end = chrono::steady_clock::now();
    // sampled->timespent = std::chrono::duration_cast<chrono::milliseconds>(end - start).count();
    // string writeFile = outputdir + "time.txt";
    // ofstream resultFile(writeFile.c_str());
    // resultFile << to_string(sampled->timespent) << endl;
    // resultFile.close();

    // Initiate
    HSet *sampled;
    set<int> initial_state;    
    for(int h = 0 ; h < graph->number_of_hedges; h++){
        htable[h] = (int)graph->hyperedge2node[h].size();
    }
    initiate();
    if (addflag){
        sampled = new HSet(initial_state, graph, outputdir, accuracy, ""); // empty
    }
    else{
        for(int h = 0 ; h < graph->number_of_hedges; h++){
            initial_state.insert(h);
        }
        sampled = new HSet(initial_state, graph, outputdir, accuracy, ""); // full
        sampled->change_version(addflag);
    }
    
    // Sample Nodes
    int iter = 0;
    while(iter < graph->number_of_hedges){
        if ((int)pool.size() > 0){
            pool.clear();
        }
        int sampled_node = sample_node();
        shuffle(graph->node2hyperedge[sampled_node].begin(), graph->node2hyperedge[sampled_node].end(), g);
        for(auto h : graph->node2hyperedge[sampled_node]){
            if (addflag){
                if (htable[h] > 0){
                    htable[h]--;
                    if(htable[h] == 0){
                        pool.push_back(h);
                    }
                }
            }
            else if (sampled->hypergraph_masking[h] == 1){
                pool.push_back(h);
            }
        }
        if ((int)pool.size() > 0){
            if (addflag){
                sampled->update(pool, graph, "+");
            }
            else{
                sampled->update(pool, graph, "-");
            }
            iter += (int)pool.size();
        }
    }
    sampled->save_as_txt(graph);
    return sampled;
}