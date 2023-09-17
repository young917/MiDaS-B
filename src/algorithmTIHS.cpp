# include "algorithmTIHS.hpp"
int Algorithm_TIHS::sample_hedge(vector<bool> hedge_check, int remain){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, remain - 1);
    int random_index = dist(gen);
    int sampled_hedge = -1;
    int idx = 0;
    for (int h = 0 ; h < (int)hedge_check.size() ; h++){
        if (!hedge_check[h]){
            if (idx == random_index){
                sampled_hedge = h;
                break;
            }
            idx += 1;
        }
    }
    return sampled_hedge;
}

HSet* Algorithm_TIHS::run(double accuracy){
    std::random_device rd;
    std::mt19937 gen(rd());
    // double thr = 1.0 / accuracy;

    string TimeOutputname = outputdir + "time_sampling_portion.txt";
    if (file_exist(TimeOutputname)){
        remove(TimeOutputname.c_str());
    }
    auto start_sampling_time = std::chrono::steady_clock::now();
    int unit_numhedge = (int)ceil((double)graph->number_of_hedges / 10.0);

    // initialize
    set<int> initial_state;
    HSet *sampled = new HSet(initial_state, graph, outputdir, accuracy, "");
    fill(hedge_check.begin(), hedge_check.end(), false);
    fill(check_node.begin(), check_node.end(), false);
    fill(check_hyperedge.begin(), check_hyperedge.end(), false);

    // start sampling
    for(int h = 0 ; h < graph->number_of_hedges ; h++){
        check_hyperedge[h] = (int)graph->hyperedge2node[h].size();
    }
    vector<int> tmp;
    while(sampled->number_of_hedges < graph->number_of_hedges){
        if ((int)pool.size() > 0){
            pool.clear();
        }
        int h = sample_hedge(hedge_check, graph->number_of_hedges - sampled->number_of_hedges);
        pool.push_back(h);
        hedge_check[h] = true;
        // cout << "sample new hyperedge" << endl;
        int hsize = (int)graph->hyperedge2node[h].size();
        for(int i = 0 ; i < hsize ; i++){
            int v = graph->hyperedge2node[h][i];
            if (!check_node[v]){
                int vdeg = (int)graph->node2hyperedge[v].size();
                shuffle(graph->node2hyperedge[v].begin(), graph->node2hyperedge[v].end(), gen);
                for (int j = 0 ; j < vdeg ; j++){
                    int h_p = graph->node2hyperedge[v][j];
                    check_hyperedge[h_p] -= 1;
                    if (check_hyperedge[h_p] == 0){
                        if (h_p != h){
                            pool.push_back(h_p);
                        }
                        hedge_check[h_p] = true;
                    }
                }
                check_node[v] = true;
            }
        }
        for (int pi=0 ; pi < (int)pool.size() ; pi++){
            tmp.push_back(pool[pi]);
            sampled->update(tmp, graph, "+");
            if ((sampled->number_of_hedges % unit_numhedge == 0) || (sampled->number_of_hedges == graph->number_of_hedges)){
                auto end_sampling_time = std::chrono::steady_clock::now();
                const auto runtime =  std::chrono::duration_cast<chrono::milliseconds>(end_sampling_time - start_sampling_time);
                ofstream TimeOutput(TimeOutputname.c_str(), std::ios::app | std::ios::out );
                TimeOutput << to_string(runtime.count()) << " ms\n";
                TimeOutput.close();
            }
            tmp.clear();
        }
    }
    sampled->save_as_txt(graph);
    
    return sampled;
}