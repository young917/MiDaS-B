#include "algorithmESSZ.hpp"

void AlgorithmESSZ::update_prop(int &hprop, int cand){
    int length = (int)algo_opt.size();
    if ( (length > 3) && (algo_opt.compare(length-3 , 3, "min") == 0) ){
        if ( (hprop == -1) || (hprop > cand) ){
            hprop = cand;
        }
    }
    else if ( (length > 3) && (algo_opt.compare(length-3 , 3, "max") == 0) ){
        if ( (hprop == -1) || (hprop < cand) ){
            hprop = cand;
        }
    }
    else if ( (length > 3) && (algo_opt.compare(length-3 , 3, "avg") == 0) ){
        hprop += cand;
    }
}

void AlgorithmESSZ::initiate(void){
    // clean data
    for (auto &item: leaf2hedges){
        auto &vec = item.second;
        vec.clear();
    }
    leaf2hedges.clear();
    hedge_tree.clear();
    max_key = 0;
    int max_hdeg = 0;
        
    // make hedge2prop by global_deg_{avg/min/max}^{alpha} * hedge_size^{-beta}
    if(algo_opt.compare(0, 10, "global_deg")==0){
        for(int v = 0 ; v < graph->number_of_nodes ; v++){
            int deg = (int)graph->node2hyperedge[v].size();
            node2prop[v] = deg;
        }
        for (int h = 0 ; h < graph->number_of_hedges ; h++){
            int hprop = -1;
            int hsize = (int)graph->hyperedge2node[h].size();
            if ( (algo_opt_length > 3) && (algo_opt.compare(algo_opt_length - 3, 3, "avg") == 0) ){
                hprop = 0;
            }
            for (int vidx = 0 ; vidx < hsize ; vidx++){
                int v = graph->hyperedge2node[h][vidx];
                update_prop(hprop, node2prop[v]);
            }
            if ( (algo_opt_length > 3) && (algo_opt.compare(algo_opt_length - 3, 3, "avg") == 0) ){
                hprop = (int)round((double)hprop / hsize);
            }
            string key = to_string(hprop) + "_" + to_string(hsize);
            int key_index;
            std::unordered_map<string, int>::iterator it = prop2index.find(key);
            if (it == prop2index.end()){
                key_index = max_key;
                prop2index[key] = key_index;
                index2prop[key_index] = key;
                max_key += 1;
            }
            else{
                key_index = it->second;
            }
            hedge2prop[h] = key_index;
            max_hdeg = max(max_hdeg, hprop);
        }
    }
    cout << "Max Key = " << to_string(max_key) << " Max Hdeg = " << to_string(max_hdeg) << endl;
    // make tree
    int height = (int)ceil(log2(max_key + 1));
    int whole_size = exp2(height+1);
    leaf_start = exp2(height);    
    hedge_tree.resize(whole_size, 0);
    for(int h = 0; h < graph->number_of_hedges ; h++){
        int hprop = hedge2prop[h];
        leaf2hedges[hprop].push_back(h);
    }

    for(int k = 0 ; k < max_key ; k++){
        if ((int)leaf2hedges[k].size() == 0){
            hedge_tree[leaf_start + k] = 0;
        }
        else{
            string hprop = index2prop[k];
            vector<string> tmp = split(hprop, '_');
            int hdeg = stoi(tmp[0]);
            int hsize = stoi(tmp[1]);
            hedge_tree[leaf_start + k] = pow(max((double)hdeg / resolution, EPSILON), alpha) * pow(max((double)hsize, 1.0), beta) * (double)leaf2hedges[k].size();
            // cout << to_string(hdeg) << " " << to_string(hsize) << " : " << to_string(hedge_tree[leaf_start + k]) << endl;
            assert (hedge_tree[leaf_start + k] > 0.0);
        }
    }
    for(int p = leaf_start - 1 ; p > 0 ; p--){
        hedge_tree[p] = hedge_tree[2 * p] + hedge_tree[2 * p + 1];
        assert (hedge_tree[p] >= 0.0);
    }
    // cout << "Start " << to_string(hedge_tree[1]) << endl;
    return;
}

int AlgorithmESSZ::sample_hedge(void){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);
    int idx = 1;
    int option_length = (int)algo_opt.size();
    while (idx < leaf_start){
        if (hedge_tree[idx] <= 0.0){
            cout << "Zero Node " << to_string(idx) << endl;
            cout << to_string(hedge_tree[(int)(round(idx / 2))]) << endl;
        }
        assert (hedge_tree[idx] > 0.0);
        if (hedge_tree[2 * idx] == 0.0){
            idx = 2 * idx + 1;
            continue;
        }
        else if (hedge_tree[2 * idx + 1] == 0.0){
            idx = 2 * idx;
            continue;
        }
        double random_double = dist(gen);
        cout << to_string(random_double) << endl;
        assert ((hedge_tree[2 * idx] + hedge_tree[2 * idx + 1]) > 0.0);
        double prob = (double)hedge_tree[2 * idx] / (hedge_tree[2 * idx] + hedge_tree[2 * idx + 1]);
        if (random_double < prob) idx = 2 * idx;
        else idx = 2 * idx + 1;
    }
    int final_index = idx - leaf_start;
    int numedges = (int)leaf2hedges[final_index].size();
    if (numedges == 0){
        cout << "Empty leaf2hedges" << endl;
        cout << to_string(hedge_tree[idx]) << " " << final_index << endl; 
        cout << to_string(hedge_tree[1]) << endl;
    }
    assert (numedges > 0);
    std::uniform_int_distribution<> dist2(0, numedges-1);
    int random_index = dist2(gen);
    int sampled_hedge = leaf2hedges[final_index][random_index];
    // erase sampled hyperedge
    // update leaf2hedges
    leaf2hedges[final_index].erase(leaf2hedges[final_index].begin() + random_index);
    if ((int)leaf2hedges[final_index].size() == 0){
        hedge_tree[leaf_start + final_index] = 0.0;
    }
    else{
        string hprop = index2prop[final_index];
        vector<string> tmp = split(hprop, '_');
        int hdeg = stoi(tmp[0]);
        int hsize = stoi(tmp[1]);
        hedge_tree[leaf_start + final_index] = pow((double)hdeg, alpha) * pow(max((double)hsize, 1.0), beta) * (double)leaf2hedges[final_index].size();
    }

    // update hedge tree
    int p = leaf_start + final_index;
    while (p > 1){
        int parent = (int)(floor(p/2));
        hedge_tree[parent] = hedge_tree[2 * parent] + hedge_tree[2 * parent + 1];
        p = parent;
        assert (hedge_tree[parent] >= 0.0);
    }
    
    return sampled_hedge;
}

HSet* AlgorithmESSZ::run(double accuracy){
    // cout << "Run ES" << endl;
    // if (addflag){
        // cout << "Add" << endl;
    // }
    // else{
        // cout << "Remove" << endl;
    // }=
    // cout << "Algorithm Option = " << algo_opt << endl;
    // cout << "Alpha = " << to_string(alpha) << endl;
    
    HSet *sampled;
    set<int> initial_state;
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
    // cout << "Done Initialization => Sampling" << endl;
    int iter = 0;
    while(iter < graph->number_of_hedges){
        int selected_hyperedge = sample_hedge();
        if ((int)pool.size() > 0){
            pool.clear();
        }
        pool.push_back(selected_hyperedge);
        if (addflag){
            sampled->update(pool, graph, "+");
        }
        else{
            sampled->update(pool, graph, "-");
        }
        iter ++;
    }
    // cout << "Done Sampling" << endl;
    sampled->save_as_txt(graph);
    return sampled;
}