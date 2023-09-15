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

    auto start = chrono::steady_clock::now();
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
    auto end = chrono::steady_clock::now();
    sampled->timespent = std::chrono::duration_cast<chrono::milliseconds>(end - start).count();
    
    while(sampled->number_of_hedges < graph->number_of_hedges){
        start = chrono::steady_clock::now();
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
                        // cout << "check h_p" << endl;
                        if (h_p != h){
                            pool.push_back(h_p);
                        }
                        // cout << "insert into pool" << endl;
                        hedge_check[h_p] = true;
                    }
                }
                check_node[v] = true;
            }
        }
        end = chrono::steady_clock::now();
        sampled->timespent += std::chrono::duration_cast<chrono::milliseconds>(end - start).count();
        // cout << "update" << endl;

        sampled->update(pool, graph, "+");
        if (sampled->number_of_nodes < 10){
            continue;
        }
        // else if (((double)sampled->number_of_nodes / graph->number_of_nodes) > thr){
        //     thr += 1.0 / accuracy;
        //     sampled->save_properties(graph);
        // }
    }
    sampled->save_as_txt(graph);
    string writeFile = outputdir + "time.txt";
    ofstream resultFile(writeFile.c_str());
    resultFile <<  to_string(sampled->timespent) << endl;
    resultFile.close();

    return sampled;
}