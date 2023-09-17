#include "algorithmES.hpp"

void AlgorithmES::update_prop(int &hprop, int cand){
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

void AlgorithmES::initiate(void){
    if (algo_opt.compare(0, 5, "decsz") == 0){
        vector<pair<int,int>> tmp;
        for (int h = 0 ; h < graph->number_of_hedges ; h++){
            int hsize = (int)graph->hyperedge2node[h].size();
            tmp.push_back(make_pair(hsize, h));
        }
        sort(tmp.begin(), tmp.end());
        for (int hi = 0 ; hi < graph->number_of_hedges ; hi++){
            int h = tmp[hi].second;
            hedges.push_back(h);
        }
    }
    else if (algo_opt.compare(0, 6, "random") == 0){
        std::random_device rd;
        std::mt19937 g(rd());
        for (int h = 0 ; h < graph->number_of_hedges ; h++){
            hedges.push_back(h);
        }
        shuffle(hedges.begin(), hedges.end(), g);
    }
    else{
        // clean data
        for (auto &item: leaf2hedges){
            auto &vec = item.second;
            vec.clear();
        }
        leaf2hedges.clear();
        hedge_tree.clear();
        // pairdegree.clear();
        max_key = -1;
        cout << "# hyperedges = " << to_string(graph->number_of_hedges) << endl;
        // make hedge2prop
        if(algo_opt.compare(0, 10, "global_deg")==0){
            for (int h = 0 ; h < graph->number_of_hedges ; h++){
                int hprop = -1;
                if ( (algo_opt_length > 3) && (algo_opt.compare(algo_opt_length - 3, 3, "avg") == 0) ){
                    hprop = 0;
                }
                int hsize = (int)graph->hyperedge2node[h].size();
                for (int vidx = 0 ; vidx < hsize ; vidx++){
                    int v = graph->hyperedge2node[h][vidx];
                    update_prop(hprop, (int)graph->node2hyperedge[v].size());
                }
                if ( (algo_opt_length > 3) && (algo_opt.compare(algo_opt_length - 3, 3, "avg") == 0) ){
                    hprop = (int)round((double)hprop / hsize);
                }
                hedge2prop[h] = hprop;
                if ( (max_key == -1) || (hprop > max_key)) {
                    max_key = hprop;
                }
            }
        }
        cout << "Max Key = " << to_string(max_key) << endl;

        // make tree
        int height = (int)ceil(log2(max_key + 1));
        int whole_size = exp2(height+1);
        leaf_start = exp2(height);    
        hedge_tree.resize(whole_size, 0);
        cout << "Make Tree" << endl;
        for(int h = 0; h < graph->number_of_hedges ; h++){
            int hprop = hedge2prop[h];
            leaf2hedges[hprop].push_back(h);
        }
        cout << "Make leaf2hedges" << endl;
        for(int k = 0 ; k <= max_key ; k++){
            if ((int)leaf2hedges[k].size() == 0){
                hedge_tree[leaf_start + k] = 0;
            }
            else{
                hedge_tree[leaf_start + k] = pow(max((double)k / resolution, EPSILON), alpha) * (double)leaf2hedges[k].size();
            }
            assert (hedge_tree[leaf_start + k] >= 0.0);
        }
        for(int p = leaf_start - 1 ; p > 0 ; p--){
            hedge_tree[p] = hedge_tree[2 * p] + hedge_tree[2 * p + 1];
            assert (hedge_tree[p] >= 0.0);
        }
        cout << "Insert values" << endl;
    }
    return;
}

int AlgorithmES::sample_hedge(void){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);
    int sampled_hedge;

    if (algo_opt.compare(0, 5, "decsz") == 0){
        sampled_hedge = hedges[index];
        index += 1;
    }
    else if (algo_opt.compare(0, 6, "random") == 0){
        sampled_hedge = hedges[index];
        index += 1;
    }
    else{
        int idx = 1;
        int option_length = (int)algo_opt.size();
        while (idx < leaf_start){
            if (hedge_tree[2 * idx] == 0){
                idx = 2 * idx + 1;
                continue;
            }
            else if (hedge_tree[2 * idx + 1] == 0){
                idx = 2 * idx;
                continue;
            }
            double random_double = dist(gen);
            assert ((hedge_tree[2 * idx] + hedge_tree[2 * idx + 1]) > 0.0);
            double prob = (double)hedge_tree[2 * idx] / (hedge_tree[2 * idx] + hedge_tree[2 * idx + 1]);
            if (random_double < prob) idx = 2 * idx;
            else idx = 2 * idx + 1;
            assert (hedge_tree[idx] > 0.0);
        }
        int final_index = idx - leaf_start;
        int numedges = (int)leaf2hedges[final_index].size();
        if (numedges == 0){
            cout << "Empty leaf2hedges" << endl;
            cout << to_string(hedge_tree[idx]) << " " << final_index << endl; 
        }
        assert (numedges > 0);
        std::uniform_int_distribution<> dist2(0, numedges-1);
        int random_index = dist2(gen);
        sampled_hedge = leaf2hedges[final_index][random_index];

        // erase sampled hyperedge
        // update leaf2hedges
        leaf2hedges[final_index].erase(leaf2hedges[final_index].begin() + random_index);
        if ((int)leaf2hedges[final_index].size() == 0){
            hedge_tree[leaf_start + final_index] = 0.0;
        }
        else if (final_index == 0){
            hedge_tree[leaf_start + final_index] = 0.0;
        }
        else{
            hedge_tree[leaf_start + final_index] = pow((double)final_index, alpha) * (double)leaf2hedges[final_index].size();
        }
        // update hedge tree
        int p = leaf_start + final_index;
        while (p > 1){
            int parent = (int)(floor(p/2));
            hedge_tree[parent] = hedge_tree[2 * parent] + hedge_tree[2 * parent + 1];
            p = parent;
            assert (hedge_tree[parent] >= 0.0);
        }
    }
    
    return sampled_hedge;
}

HSet* AlgorithmES::run(double accuracy){
    cout << "Run ES" << endl;
    if (addflag){
        cout << "Add" << endl;
    }
    else{
        cout << "Remove" << endl;
    }
    cout << "Algorithm Option = " << algo_opt << endl;
    cout << "Alpha = " << to_string(alpha) << endl;
    
    string TimeOutputname = outputdir + "time_sampling_portion.txt";
    if (file_exist(TimeOutputname)){
        remove(TimeOutputname.c_str());
    }
    auto start_sampling_time = std::chrono::steady_clock::now();
    int unit_numhedge = (int)ceil((double)graph->number_of_hedges / 10.0);

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
    cout << "Done Initialization => Sampling" << endl;
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
        if ((iter % unit_numhedge == 0) || (iter == graph->number_of_hedges)){
            auto end_sampling_time = std::chrono::steady_clock::now();
            const auto runtime =  std::chrono::duration_cast<chrono::milliseconds>(end_sampling_time - start_sampling_time);
            ofstream TimeOutput(TimeOutputname.c_str(), std::ios::app | std::ios::out);
            TimeOutput << to_string(runtime.count()) << " ms\n";
            TimeOutput.close();
        }
    }
    cout << "Done Sampling" << endl;
    sampled->save_as_txt(graph);
    return sampled;
}