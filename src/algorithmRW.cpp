# include "algorithmRW.hpp"

Algorithm_RW::Algorithm_RW(string algo_opt, string outputdir, HyperGraph *graph, double restart, double maxlength, bool noinduce){
    this->algo_opt = algo_opt;
    this->outputdir = outputdir;
    this->graph = graph;
    this->restart = restart;
    this->givenmaxlength = maxlength;
    this->noinduce = noinduce;
    this->node2node.resize(graph->number_of_nodes);
    this->check.resize(graph->number_of_nodes, false);
    this->htable.resize(graph->number_of_hedges);

    for(int va = 0 ; va < graph->number_of_nodes ; va++){
        for(auto h : graph->node2hyperedge[va]){
            int hsize = (int)graph->hyperedge2node[h].size();
            for(int i = 0 ;  i < hsize; i++){
                int vb = graph->hyperedge2node[h][i];
                if(va != vb){
                    this->node2node[va].push_back(vb);
                    this->node2node[vb].push_back(va);
                }
            }
        }
    }
}

void Algorithm_RW::walk(int seed, int max_length){
    default_random_engine e;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);
    
    q.push(seed);
    int step = 0;
    while ((step < max_length) && (q.empty() == false)){
        int current = q.front();
        q.pop();
        if(!check[current]){
            // add induced hyperedges randomly
            check[current] = true;
            int deg = (int)graph->node2hyperedge[current].size();
            if ((int)order.size() > 0){
                order.clear();
            }
            for (int hidx = 0 ; hidx < deg ; hidx ++){
                order.push_back(hidx);
            }
            shuffle(order.begin(), order.end(), gen);
            for (int i = 0 ; i < deg ; i ++){
                int hidx = order[i];
                int h = graph->node2hyperedge[current][hidx];
                if (htable[h] > 0){
                    htable[h]--;
                    if (htable[h] == 0){
                        pool.push_back(h);
                        if (noinduce){
                            // Add only one hyperedge?
                            break;
                        }
                    }
                }
            }
        }
        // Add only one hyperedge?
        if ((noinduce) && ((int)pool.size() > 0)){
            break;
        }
        double random_double = dist(gen);
        if(random_double < restart){
            q.push(seed);
        }
        else if ((int)node2node[current].size() != 0){
            int deg = (int)node2node[current].size();
            std::uniform_int_distribution<> dist2(0, deg-1);
            int random_index = dist2(gen);
            int next_node = node2node[current][random_index];
            q.push(next_node);
        }
        step++;
    }
    while(q.empty() == false){
        q.pop();
    }
}

HSet* Algorithm_RW::run(double accuracy){
    std::random_device rd;
    std::mt19937 gen(rd());

    int max_length;
    // double thr = 1.0 / accuracy;
    if (givenmaxlength == -1){
        max_length = 100 * graph->number_of_nodes;            
    }
    else if (givenmaxlength == 1){
        max_length = graph->number_of_nodes;
    }
    else{
        max_length = givenmaxlength * graph->number_of_nodes;
    }

    string TimeOutputname = outputdir + "time_sampling_portion.txt";
    if (file_exist(TimeOutputname)){
        remove(TimeOutputname.c_str());
    }
    auto start_sampling_time = std::chrono::steady_clock::now();
    int unit_numhedge = (int)ceil((double)graph->number_of_hedges / 10.0);
    
    set<int> initial_state;
    HSet *sampled = new HSet(initial_state, graph, outputdir, accuracy, "");
    auto start = chrono::steady_clock::now();
    fill(check.begin(), check.end(), false);
    htable.resize(graph->number_of_hedges);
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        htable[h] = (int)graph->hyperedge2node[h].size();
    }
    
    vector<int> tmp2;
    while(sampled->number_of_hedges < graph->number_of_hedges){
        cout << "\r" << to_string(sampled->number_of_hedges) << "/" << to_string(graph->number_of_hedges);
        start = chrono::steady_clock::now();
        int seed;
        tmp.clear();
        for (int v = 0 ; v < graph->number_of_nodes ; v++){
            if (sampled->node_masking[v] == 0){
                tmp.push_back(v);
            }
        }
        std::uniform_int_distribution<> dist(0, (int)tmp.size() - 1);
        int seed_idx = dist(gen);
        seed = tmp[seed_idx];

        if ((int)pool.size() > 0){
            pool.clear();
            if (noinduce){
                for (int h = 0 ; h < graph->number_of_hedges ; h++){
                    if (htable[h] != 0){
                        htable[h] = (int)graph->hyperedge2node[h].size();
                    }
                }
            }
        }
        walk(seed, max_length);
        if ((int)pool.size() > 0){
            // sampled->update(pool, graph, "+");
            for (int pi=0 ; pi < (int)pool.size() ; pi++){
                tmp2.push_back(pool[pi]);
                sampled->update(tmp2, graph, "+");
                if ((sampled->number_of_hedges % unit_numhedge == 0) || (sampled->number_of_hedges == graph->number_of_hedges)){
                    auto end_sampling_time = std::chrono::steady_clock::now();
                    const auto runtime =  std::chrono::duration_cast<chrono::milliseconds>(end_sampling_time - start_sampling_time);
                    ofstream TimeOutput(TimeOutputname.c_str(), std::ios::app | std::ios::out);
                    TimeOutput << to_string(runtime.count()) << " ms\n";
                    TimeOutput.close();
                }
                tmp2.clear();
            }
        }
    }
    cout << endl;
    sampled->save_as_txt(graph);

    return sampled;
}