#include "helper.hpp"

Helper::Helper(set<int> subhypergraph, HyperGraph *graph, string outputdir, string algo_opt) {
    this->graph = graph;
    this->hypergraph_masking.resize(graph->number_of_hedges, false);
    this->node_masking.resize(graph->number_of_nodes, false);
    this->nodes.resize(graph->number_of_nodes, 0);
    this->ancestor.resize(graph->number_of_nodes);
    this->group_size.resize(graph->number_of_nodes, 1);
    for (int v = 0 ; v < graph->number_of_nodes ; v++){
        this->ancestor[v] = v;
    }
    number_of_nodes = 0;
    sum_of_hsizes=0;
    // assert ((int)subhypergraph.size() == 0);
    this->number_of_hedges = subhypergraph.size();
    this->number_of_nodes = number_of_nodes;
    this->outputdir = outputdir;

    this->check.resize(graph->number_of_hedges, 0);
    this->check_node.resize(graph->number_of_nodes, false);
    // cout << "Number of Nodes = " << graph->number_of_nodes << endl;
    
    // setting flag
    if ((int)algo_opt.size() == 0){
        ist_flag = true;
        ist_size_flag = true;
        dsf_flag = true;
        cc_flag = true;
        wcc_flag = true;
    }
    else{
        tmp = split(algo_opt, ',');
        for (int i = 0 ; i < (int)tmp.size() ; i++){
            if (tmp[i].compare("intersection") == 0) ist_flag = true;
            else if (tmp[i].compare("intersect_avg") == 0) ist_size_flag = true;
            else if (tmp[i].compare("densification") == 0) dsf_flag = true;
            else if (tmp[i].compare("clusteringcoef") == 0) cc_flag = true;
            else if (tmp[i].compare("sizewcc") == 0) wcc_flag = true;
        }
    }

    if (cc_flag){
        tmp = split(outputdir, '/');
        for (int i = 0 ; i < (int)tmp.size() ; i++){
            if (((int)tmp[i].size() > 7) && (tmp[i].compare(0, 7, "threads") == 0)){
                cc_approx = true;
                dataname = "threads";
                break;
            }
            else if(((int)tmp[i].size() > 4) && (tmp[i].compare(0, 4, "tags") == 0)){
                cc_approx = true;
                dataname = "tags";
                break;
            }
            else if (((int)tmp[i].size() > 6) && (tmp[i].compare(0, 6, "coauth") == 0)){
                cc_approx = true;
                dataname = "coauth";
                break;
            }
        }
    }
    if (cc_approx){
        cout << "Calculate CC APPROX" << endl;
    }

    // update by subhypergraph
    intersect = 0;
    if ((int)subhypergraph.size() > 0){
        set<int> nodeset;
        for (int h : subhypergraph){
            fill(check.begin(), check.end(), 0);
            this->hypergraph_masking[h] = true;
            int hsize = (int)graph->hyperedge2node[h].size();
            this->sum_of_hsizes += hsize;
            for (int i = 0; i < hsize ; i++){
                int vi = graph->hyperedge2node[h][i];
                nodeset.insert(vi);
                if (!this->node_masking[vi]){
                    this->number_of_nodes++;
                    this->node_masking[vi] = true; // update node_masking
                }
                if (cc_flag || wcc_flag){
                    for (int j = i+1 ; j < hsize ; j++){
                        int vj = graph->hyperedge2node[h][j];
                        if (wcc_flag){
                            // update ancestor
                            int anc_vi = get_ancestor(vi);
                            int anc_vj = get_ancestor(vj);
                            if (anc_vi != anc_vj){
                                if (anc_vi < anc_vj){
                                    ancestor[anc_vj] = anc_vi;
                                    group_size[anc_vi] += group_size[anc_vj];
                                    group_size[anc_vj] = 0;
                                }
                                else{
                                    ancestor[anc_vi] = anc_vj;
                                    group_size[anc_vj] += group_size[anc_vi];
                                    group_size[anc_vi] = 0;
                                }
                            }
                        }
                        if (cc_flag){
                            string key = make_sortedkey(vi, vj);
                            this->neighbor.insert(key);
                        }
                    }
                    
                }
                if (ist_flag || ist_size_flag){
                    int vdeg = (int)graph->node2hyperedge[vi].size();
                    for (int j = 0 ; j < vdeg ; j++){
                        int nh = graph->node2hyperedge[vi][j];
                        if (!hypergraph_masking[nh]){
                            continue;
                        }
                        if (h != nh){
                            if ( (ist_flag) && (check[nh] == 0) ){
                                intersect += 1;
                            }
                            int prev = check[nh];
                            check[nh] += 1;
                            if ((prev > 0) && (ist_size_flag)){
                                assert (intersection_size[prev] > 0);
                                intersection_size[prev]--;
                                int cur = check[nh];
                                intersection_size[cur]++;
                            }
                        }
                        // if ((h != nh) && (check[nh] == false)){
                        //     // update intersect connection
                        //     intersect += 1;
                        //     // update intersect size
                        //     if (ist_size_flag){
                        //         vec.resize(graph->number_of_nodes);
                        //         it = set_intersection(graph->hyperedge2node[h].begin(), graph->hyperedge2node[h].end(), graph->hyperedge2node[nh].begin(), graph->hyperedge2node[nh].end(), vec.begin());
                        //         int its_size = (int)(it - vec.begin());
                        //         intersection_size[its_size] += 1;
                        //         vec.clear();
                        //     }
                        // }
                    }
                }
            }
            // if (ist_size_flag){
            //     for (int nh = 0 ; nh < graph->number_of_hedges ; nh++){
            //         if (check[nh] > 0){
            //             int its_size = check[nh];
            //             intersection_size[its_size] += 1;
            //         }
            //     }
            // }
        }
        this->number_of_nodes = (int)nodeset.size();
    }
}

void Helper::get_intersection(void){
    long long denom = (long long)number_of_hedges * (long long)(number_of_hedges - 1) / 2;
    long long numer = intersect;
    // cout << numer << " " << denom <<  endl;
    string writeFile = outputdir + "intersection.txt";
    ofstream resultFile(writeFile.c_str(), fstream::out | fstream::app);
    // resultFile << "num nodes,num edges"  << endl;
    resultFile << to_string(numer) << "," << to_string(denom) << endl;
    resultFile.close();
}

void Helper::get_intersection_size(void){
    double avgintersectsize = 0;
    long long total = 0;
    for (auto elem : intersection_size){
        total += elem.second;
    }
    for (auto elem : intersection_size){
        avgintersectsize += (double)elem.first * ((double)elem.second / (double)total);
    }
    string writeFile = outputdir + "intersect_avg.txt";
    ofstream resultFile(writeFile.c_str(), fstream::out | fstream::app);
    resultFile << to_string(number_of_hedges) << "," << to_string(avgintersectsize) << endl;
    resultFile.close();
}

void Helper::get_densification(void){
    string writeFile = outputdir + "densification.txt";
    //cout << "Output Density of graph " << outputdir << endl;
    ofstream resultFile(writeFile.c_str(), fstream::out | fstream::app);
    // resultFile << "num nodes,num edges"  << endl;
    resultFile << to_string(number_of_nodes) << "," << to_string(number_of_hedges) << endl;
    resultFile.close();
}

void Helper::get_average_clustering_coef_appx(void){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, number_of_nodes-1);
    
    double epsilon;
    double delta;
    double k;
    // probability (1- \delta) ; |x-\bar{x}| < \epsilon
    if (dataname.compare("tags") == 0){
        epsilon = 0.1;
        delta = 0.1;
        k = 0.5 * pow(epsilon, -2) * log(2.0 / delta);
    }
    else if (dataname.compare("coauth") == 0){
        epsilon = 0.1;
        delta = 0.1;
        k = 0.5 * pow(epsilon, -2) * log(2.0 / delta);
        // cout << "K is " << to_string(k) << endl;
    }
    else if (dataname.compare("threads") == 0){
        epsilon = 0.1;
        delta = 0.1;
        k = 0.5 * pow(epsilon, -2) * log(2.0 / delta);
        // cout << "K is " <<  to_string(k) << endl;
    }

    int step = 0;
    int closed = 0;
    string key;
    vector<int> allnodes;
    for (int v = 0 ; v < graph->number_of_nodes; v++){
        if (node_masking[v]) allnodes.push_back(v);
    }
    while (step < k){
        step += 1;
        int random_index = dist(gen);
        int v = allnodes[random_index];
        int i = 0;
        for (int nv = 0 ; nv < (int)node_masking.size() ; nv++){
            if ((!node_masking[nv]) || (v == nv)){
                continue;
            }
            key = make_sortedkey(v, nv);
            if (neighbor.find(key) != neighbor.end()){
                nodes[i] = nv;
                i++;
            }
        }
        int vdeg = i;
        if (vdeg >= 2){
            shuffle(nodes.begin(), nodes.begin() + vdeg, std::default_random_engine());
            int va = nodes[0];
            int vb = nodes[1];
            key = make_sortedkey(va, vb);
            if (neighbor.find(key) != neighbor.end()){
                closed++;
            }
        }
    }
    double average = (double)closed / (double)step;
    string writeFile = outputdir + "clusteringcoef.txt";
    ofstream resultFile(writeFile.c_str(), fstream::out | fstream::app);
    resultFile << to_string(number_of_hedges) << "," << to_string(average) << "," <<  to_string(number_of_nodes) << endl;
    resultFile.close();
}

void Helper::get_average_clustering_coef(void){
    double average = 0;
    string key;

    for (int v = 0; v < (int)node_masking.size() ; v++){
        if (!node_masking[v]){
            continue;
        }
        int i = 0;
        for (int nv = 0 ; nv < (int)node_masking.size() ; nv++){
            if ((!node_masking[nv]) || (v == nv)){
                continue;
            }
            key = make_sortedkey(v, nv);
            if (neighbor.find(key) != neighbor.end()){
                nodes[i] = nv;
                i++;
            }
        }
        int vdeg = i;
        double cc = 0.0; // number of connected neighbor pairs
        double denominator = 0.0; // number of neighbor pairs
        if (vdeg < 2){
            continue;
        }
        // neighbor pair
        for(int va_idx = 0 ; va_idx < vdeg ; va_idx++){
            for(int vb_idx = va_idx + 1 ; vb_idx < vdeg ; vb_idx++){
                int va = nodes[va_idx];
                int vb = nodes[vb_idx];
                key = make_sortedkey(va, vb);
                if (neighbor.find(key) != neighbor.end()){
                    cc++;
                }
                denominator++;
            }
        }
        average += (cc / denominator) / number_of_nodes;
        // cout << to_string(cc) << " " << to_string(denominator) << " " << to_string(average) << endl;
    }
    string writeFile = outputdir + "clusteringcoef.txt";
    ofstream resultFile(writeFile.c_str(), fstream::out | fstream::app);
    resultFile << to_string(number_of_hedges) << "," << to_string(average) << "," <<  to_string(number_of_nodes) << endl;
    resultFile.close();
}

int Helper::get_ancestor(int node){
    if (ancestor[node] == node) return node;
    else return get_ancestor(ancestor[node]);
}

int Helper::bfs(int start_node){
    int number_of_visited = 0;
    visited.push(start_node);
    check_node[start_node] = true;
    number_of_visited++;
    while((int)visited.size() > 0){
        int v = visited.front();
        visited.pop();
        for (int nv = 0 ; nv < (int)node_masking.size() ; nv++){
            if ((!node_masking[nv]) || (v == nv)){
                continue;
            }
            string key = make_sortedkey(v, nv);
            if (neighbor.find(key) != neighbor.end()){
                if (check_node[nv] == false){
                    visited.push(nv);
                    check_node[nv] = true;
                    number_of_visited++;
                }
            }
        }
    }
    return number_of_visited;
}

void Helper::count_wcc(void){
    // int largest_wcc_size = 0;
    // fill(check_node.begin(), check_node.end(), false);
    // for(int start = 0 ; start < (int)node_masking.size() ; start++){
    //     if( (node_masking[start]) && (check_node[start] == false)){
    //         int size_wcc = bfs(start);
    //         if (largest_wcc_size < size_wcc){
    //             largest_wcc_size = size_wcc;
    //         }
    //     }
    // }
    // string writeFile = outputdir + "sizewcc.txt";
    // ofstream resultFile(writeFile.c_str(), fstream::out | fstream::app);
    // resultFile << to_string(number_of_hedges) << "," << to_string(largest_wcc_size) << "," << to_string(number_of_nodes) << endl;
    // resultFile.close();

    int largest_wcc_size = 0;
    for (int v = 0 ; v < graph->number_of_nodes; v++){
        if (group_size[v] > largest_wcc_size) largest_wcc_size = group_size[v];
    }
    // assert (new_largest_wcc_size == largest_wcc_size);
    string writeFile = outputdir + "sizewcc.txt";
    ofstream resultFile(writeFile.c_str(), fstream::out | fstream::app);
    resultFile << to_string(number_of_hedges) << "," << to_string(largest_wcc_size) << "," << to_string(number_of_nodes) << endl;
    resultFile.close();
}

void Helper::get_dense_property(void){
    double density = double(number_of_hedges) / number_of_nodes;
    double overlapness = sum_of_hsizes / number_of_nodes;

    string writeFile1 = outputdir + "density.txt";
    ofstream resultFile1(writeFile1.c_str(),  fstream::out | fstream::app);
    resultFile1 << to_string(number_of_hedges) << "," << to_string(density) << endl;
    resultFile1.close();

    string writeFile2 = outputdir + "overlapness.txt";
    ofstream resultFile2(writeFile2.c_str(),  fstream::out | fstream::app);
    resultFile2 << to_string(number_of_hedges) << "," << to_string(overlapness) << endl;
    resultFile2.close();
}

void Helper::update(set<int> deltaset, HyperGraph *graph){
    for (auto h : deltaset){
        assert(!hypergraph_masking[h]);
        hypergraph_masking[h] = true;
        int hsize = (int)graph->hyperedge2node[h].size();
        sum_of_hsizes += hsize;
        number_of_hedges++;

        for (int i = 0 ; i < hsize ; i++){
            int v = graph->hyperedge2node[h][i];
            if (!node_masking[v]){
                number_of_nodes++;
                node_masking[v] = true;
            }
        }
        // intersect
        if (ist_flag || ist_size_flag){
            fill(check.begin(), check.end(), false);
            for (int i = 0 ; i < hsize ; i++){
                int v = graph->hyperedge2node[h][i];
                int vdeg = (int)graph->node2hyperedge[v].size();
                for (int j = 0 ; j < vdeg ; j++){
                    int nh = graph->node2hyperedge[v][j]; // neighbor hyperedge
                    if (!hypergraph_masking[nh]){
                        continue;
                    }
                    if (h != nh){
                        if ( (ist_flag) && (check[nh] == 0) ){
                            intersect += 1;
                        }
                        int prev = check[nh];
                        check[nh] += 1;
                        if (ist_size_flag){
                            if (prev > 0){
                                assert (intersection_size[prev] > 0);
                                intersection_size[prev]--;
                            }
                            int cur = check[nh];
                            intersection_size[cur]++;
                        }
                    }
                    // if ((h != nh) && (check[nh] == false)){
                    //     check[nh] = true;
                    //     // 1. connection
                    //     intersect += 1;
                    //     // 2. intersection size
                    //     if (ist_size_flag){
                    //         vec.resize(graph->number_of_nodes);
                    //         it = set_intersection(graph->hyperedge2node[h].begin(), graph->hyperedge2node[h].end(), graph->hyperedge2node[nh].begin(), graph->hyperedge2node[nh].end(), vec.begin());
                    //         int its_size = (int)(it - vec.begin());
                    //         intersection_size[its_size] += 1;
                    //         vec.clear();
                    //     }
                    // }
                    assert(intersect >= 0);
                }
            }
            // if (ist_size_flag){
            //     for (int nh = 0 ; nh < graph->number_of_hedges ; nh++){
            //         if (check[nh] > 0){
            //             int its_size = check[nh];
            //             intersection_size[its_size] += 1;
            //         }
            //     }
            // }
        }
        // neighbors
        if (cc_flag || wcc_flag){
            for (int i = 0 ; i < hsize ; i++){
                for ( int j = i+1 ; j < hsize ; j++){
                    int v = graph->hyperedge2node[h][i];
                    int nv = graph->hyperedge2node[h][j];
                    if (wcc_flag){
                        int anc_v = get_ancestor(v);
                        int anc_nv = get_ancestor(nv);
                        if (anc_v != anc_nv){
                            if (anc_v < anc_nv){
                                ancestor[anc_nv] = anc_v;
                                group_size[anc_v] += group_size[anc_nv];
                                group_size[anc_nv] = 0;
                            }
                            else{
                                ancestor[anc_v] = anc_nv;
                                group_size[anc_nv] += group_size[anc_v];
                                group_size[anc_v] = 0;
                            }
                        }
                    }
                    if (cc_flag){
                        string key = make_sortedkey(v, nv);
                        neighbor.insert(key);
                    }
                }
            }
        }
    }
    return;
}

void Helper::save_properties(void){
    if (ist_flag) get_intersection();
    if (ist_size_flag) get_intersection_size();
    if (dsf_flag) get_densification();
    if (cc_flag){
        if(cc_approx && ((int)neighbor.size() > 5000)) get_average_clustering_coef_appx(); 
        else get_average_clustering_coef();
    }
    if (wcc_flag) count_wcc();
}