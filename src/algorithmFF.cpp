# include "algorithmFF.hpp"

Algorithm_FF::Algorithm_FF(double p, double q, string algo_opt, string outputdir, HyperGraph *graph){
    this->p = p;
    this->q = q;
    this->algo_opt = algo_opt;
    this->outputdir = outputdir;
    this->graph = graph;
    this->htable.resize(graph->number_of_hedges, 0);
    this->check_node.resize(graph->number_of_nodes, false);
    // this->tie.resize(graph->number_of_nodes, vector<int>(graph->number_of_nodes, 0));
    // for (int h = 0; h < graph->number_of_hedges ; h++){
    //     int hsize = graph->hyperedge2node[h].size();
    //     for (int i = 0; i < hsize; i++){
    //         for(int j = i + 1 ; j < hsize; j++){
    //             int va = graph->hyperedge2node[h][i];
    //             int vb = graph->hyperedge2node[h][j];
    //             tie[va][vb] += 1;
    //         }
    //     }
    // }
}

void Algorithm_FF::burn1(int ambassador, vector<int> &burned_nodes_list, double prob){
    random_device rd;
    mt19937 gen(rd());
    geometric_distribution<> d(1 - prob);
    default_random_engine e;
    vector<bool> burned_nodes(graph->number_of_nodes);
    vector<bool> inqueue(graph->number_of_nodes);
    vector<pair<int, int>> candidate_nodes;
    vector<int> tie(graph->number_of_nodes);

    if(algo_opt.compare("ff_c") == 0){
        queue<int> q;
        q.push(ambassador);
        inqueue[ambassador] = true;
        bool exit = false;
        while (!q.empty()){
            int v = q.front();
            q.pop();
            burned_nodes[v] = true;
            burned_nodes_list.push_back(v);
            inqueue[v] = false;

            fill(tie.begin(), tie.end(), 0);
            int vdeg = (int)graph->node2hyperedge[v].size();
            for (int hidx = 0 ; hidx < vdeg ; hidx++){
                int h = graph->node2hyperedge[v][hidx];
                int hsize = (int)graph->hyperedge2node[h].size();
                for (int nvidx = 0; nvidx < hsize ; nvidx++){
                    int nv = graph->hyperedge2node[h][nvidx];
                    if (v != nv){
                        tie[nv]++;
                    }
                }
            }
            int n = d(gen);
            candidate_nodes.clear(); // fill with unvisited neighbor nodes
            for(int nv = 0 ; nv < graph->number_of_nodes ; nv++){
                if (((burned_nodes[nv] == false) && (inqueue[nv] == false))){
                    candidate_nodes.push_back(make_pair(-tie[nv], nv));
                }
            }
            sort(candidate_nodes.begin(), candidate_nodes.end());
            for( int i = 0 ; i < min(n, (int)candidate_nodes.size()) ; i++){ // visit n ~ geometric dist. neighbors
                int next_node = candidate_nodes[i].second;
                q.push(next_node);
                inqueue[next_node] = true;
            }
        }
    }
}

void Algorithm_FF::burn2(int ambassador, set<int> &add, double prob){
    random_device rd;
    mt19937 gen(rd());
    geometric_distribution<> d(1 - prob);
    default_random_engine e;
    vector<bool> burned_nodes(graph->number_of_nodes);
    vector<bool> inqueue(graph->number_of_nodes);
    vector<pair<int, int>> candidate_nodes;
    vector<int> tie(graph->number_of_nodes);

    if(algo_opt.compare("ff_c") == 0){
        queue<int> q;
        q.push(ambassador);
        inqueue[ambassador] = true;
        bool exit = false;
        while (!q.empty()){
            int v = q.front();
            q.pop();
            burned_nodes[v] = true;
            if (check_node[v] == false){
                check_node[v] = true;
                int deg = graph->node2hyperedge[v].size();
                tmp.clear();
                for (int hidx = 0 ; hidx < deg ; hidx ++){
                    tmp.push_back(hidx);
                }
                shuffle(tmp.begin(), tmp.end(), gen);
                for (int i = 0 ; i < deg ; i++){
                    int hidx = tmp[i];
                    int h = graph->node2hyperedge[v][hidx];
                    htable[h]--;
                    if (htable[h] == 0){
                        add.insert(h);
                        cout << "\r" << to_string((int)add.size());
                    }
                }
            }
            inqueue[v] = false;
            fill(tie.begin(), tie.end(), 0);

            int vdeg = (int)graph->node2hyperedge[v].size();
            for (int hidx = 0 ; hidx < vdeg ; hidx++){
                int h = graph->node2hyperedge[v][hidx];
                int hsize = (int)graph->hyperedge2node[h].size();
                for (int nvidx = 0; nvidx < hsize ; nvidx++){
                    int nv = graph->hyperedge2node[h][nvidx];
                    if (v != nv){
                        tie[nv]++;
                    }
                }
            }
            int n = d(gen);
            candidate_nodes.clear(); // fill with unvisited neighbor nodes
            for(int nv = 0 ; nv < graph->number_of_nodes ; nv++){
                if (((burned_nodes[nv] == false) && (inqueue[nv] == false))){
                    candidate_nodes.push_back(make_pair(-tie[nv], nv));
                }
            }
            sort(candidate_nodes.begin(), candidate_nodes.end());
            for( int i = 0 ; i < min(n, (int)candidate_nodes.size()) ; i++){ // visit n ~ geometric dist. neighbors
                int next_node = candidate_nodes[i].second;
                q.push(next_node);
                inqueue[next_node] = true;
            }
        }
    }
}

HSet* Algorithm_FF::run(double accuracy){
    std::random_device rd;
    std::mt19937 gen(rd());
    // std::uniform_int_distribution<> dist(0, graph->number_of_nodes - 1);

    string TimeOutputname = outputdir + "time_sampling_portion.txt";
    if (file_exist(TimeOutputname)){
        remove(TimeOutputname.c_str());
    }
    auto start_sampling_time = std::chrono::steady_clock::now();
    int unit_numhedge = (int)ceil((double)graph->number_of_hedges / 10.0);

    set<int> initial_state;
    HSet *sampled = new HSet(initial_state, graph, outputdir, accuracy, "");
    // initialize htalbe and check_node
    for (int h = 0 ; h < graph->number_of_hedges ; h++){
        htable[h] = (int)graph->hyperedge2node[h].size();
    }
    fill(check_node.begin(), check_node.end(), false);
    
    set<int> add;
    vector<int> pool;
    while(sampled->number_of_hedges < graph->number_of_hedges){
        tmp.clear();
        for (int v = 0 ; v < graph->number_of_nodes ; v++){
            if (sampled->node_masking[v] == 0){
                tmp.push_back(v);
            }
        }
        std::uniform_int_distribution<> dist(0, (int)tmp.size() - 1);
        int seed_idx = dist(gen);
        int ambassador = tmp[seed_idx];
        // int ambassador = dist(gen);
        
        tmp.clear();
        pool.clear();
        burn1(ambassador, pool, p);
        for(int i = 0 ; i < (int)pool.size() ; i++){
            int next_ambassador = pool[i];
            burn2(next_ambassador, add, q);
            cout << endl;
            cout << "update sample" << endl;
            // sampled->update(add, graph, "+");
            tmp.clear();
            for (auto h : add){
                tmp.push_back(h);
                sampled->update(tmp, graph, "+");
                if ((sampled->number_of_hedges % unit_numhedge == 0) || (sampled->number_of_hedges == graph->number_of_hedges)){
                    cout << sampled->number_of_hedges << endl;
                    auto end_sampling_time = std::chrono::steady_clock::now();
                    const auto runtime =  std::chrono::duration_cast<chrono::milliseconds>(end_sampling_time - start_sampling_time);
                    ofstream TimeOutput(TimeOutputname.c_str(), std::ios::app | std::ios::out);
                    TimeOutput << to_string(runtime.count()) << " ms\n";
                    TimeOutput.close();
                }
                tmp.clear();
            }
            cout << sampled->number_of_hedges << " / " << graph->number_of_hedges << endl;
            add.clear();
        }
    }
    cout << endl;
    sampled->save_as_txt(graph);
    
    return sampled;
}