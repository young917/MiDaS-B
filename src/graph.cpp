#include "graph.hpp"

HyperGraph::HyperGraph(string inputpath, string dataname){
    string path;
    this->inputpath = inputpath;
    this->dataname = dataname;
    if( (dataname.compare("cora") == 0) || (dataname.compare("citeseer") == 0)){
        this->exist_edgename = true;
        path = inputpath + dataname + "/" + dataname + ".hyperedge";
        cout << "Datapath is " << path << endl;
    }
    else if (inputpath.compare(0, 7, "results") == 0){
        path = inputpath + dataname + "/2/sampled_graph.txt";
        this->exist_edgename = false;
    }
    else{
        cout << path << endl;
        path = inputpath + dataname + "_final.txt";
        this->exist_edgename = false;
    }

	ifstream graphFile(path.c_str());
	string line;
	int num_hyperedge = 0;
    int num_nodes = 0;
    cout << "Start Read Dataset " << path << endl; 
	while (getline(graphFile, line)){
        string ename;
        if (this->exist_edgename){
            vector<string> tmp = split(line, '\t');
            ename = tmp[0];
            line = tmp[1];
        }
		vector<string> nodes = split(line, ',');
		vector<int> tokens;
		for (int i = 0; i < (int)nodes.size(); i++){
			// tokens.push_back(stoi(nodes[i]));
            if (nodename2index.find(nodes[i]) == nodename2index.end()){
                int index = num_nodes++;
                nodename2index[nodes[i]] = index;
                index2nodename[index] = nodes[i];
                this->node2hyperedge.push_back(vector<int>());
            }
            int node_index = nodename2index[nodes[i]];
            tokens.push_back(node_index);
            this->node2hyperedge[node_index].push_back(num_hyperedge);
        }
        sort(tokens.begin(), tokens.end());
        this->hyperedge2node.push_back(tokens);
        if (this->exist_edgename){
            edgename[num_hyperedge] = ename;
        }
        num_hyperedge++;
	}
    this->number_of_hedges = num_hyperedge;
    this->number_of_nodes = (int)this->node2hyperedge.size();
    cout << "Load " << number_of_hedges << " hyperedges" << endl;

    index2order.resize(number_of_nodes, -1);
    for (int h =  0 ; h < number_of_hedges ; h++){
        int hsize = hyperedge2node[h].size();
        for (int vi= 0 ; vi < hsize; vi++){
            int v = hyperedge2node[h][vi];
            if (index2order[v] == -1){
                index2order[v] = h;
            }
        }
    }
}
vector<vector<int>> HyperGraph::get_incidence_matrix(){
    int row = this->number_of_hedges;
    int col = this->number_of_nodes;
    vector<vector<int>> inc_mat(row, std::vector<int>(col, 0));
    for (int h = 0 ; h < row ; h++){
        for (auto v : this->hyperedge2node[h]){
            inc_mat[h][v] = 1;
        }
    }
    return inc_mat;
}