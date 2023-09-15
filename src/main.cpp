#include <sys/types.h>
#include <ctime>

#include "algorithmNS.hpp"
#include "algorithmES.hpp"
#include "algorithmESSZ.hpp"
#include "algorithmFF.hpp"
#include "algorithmRW.hpp"
#include "algorithmTIHS.hpp"
#include "helper.hpp"
#include "helperdist.hpp"

using namespace std;

int main(int argc, char* argv[]){
    string inputpath;
    string dataname;
    // algorihtm options
    string algorithm;
    string algo_opt = "";
    string inputdir = "";
    double accuracy = 500;
    double samplingportion = 0.5;
    // ES
    double alpha = 1.0;
    double beta = 1.0;
    // RW or FF options
    double maxlength = -1;
    double restart = 0.15;
    double p = 0.51;
    double q = 0.20;
    bool noinduce = false;
    // Flags
    bool recalculate = false;
    // Test
    int repeat = -1;

    for(int i=1; i<argc ; i++){
        string input = argv[i];
        if(input.compare("--inputpath") == 0) inputpath = argv[++i];
        else if(input.compare("--dataname") == 0) dataname = argv[++i];
        else if(input.compare("--accuracy") == 0) accuracy = atof(argv[++i]);
        else if(input.compare("--samplingportion") == 0) samplingportion = atof(argv[++i]);
        else if(input.compare("--inputdir") == 0) inputdir = argv[++i];

        else if(input.compare("--algorithm") == 0) algorithm = argv[++i];
        else if(input.compare("--algo_opt") == 0) algo_opt = argv[++i]; 
        // es
        else if (input.compare("--alpha") == 0) alpha = atof(argv[++i]);
        else if (input.compare("--beta") == 0) beta = atof(argv[++i]);
        // rw
        else if(input.compare("--noinduce") == 0) noinduce = true;
        else if(input.compare("--maxlength") == 0) maxlength = atof(argv[++i]);
        else if(input.compare("--restart") == 0) restart = atof(argv[++i]);
        // ff & ffs
        else if(input.compare("--p") == 0) p = atof(argv[++i]);
        else if(input.compare("--q") == 0) q = atof(argv[++i]);
        // run option
        else if(input.compare("--recalculate") == 0) recalculate = true;
        // test for dynamic computation
        else if (input.compare("--repeat") == 0) repeat = atoi(argv[++i]);
    }
    HyperGraph *graph = new HyperGraph(inputpath, dataname);
    string outputdir = "";
    // string portion_str = "";
    vector<string> args;
    args.push_back(algorithm);
    if((algorithm.compare(0, 6, "answer") != 0) and ((int)algo_opt.size() != 0)){
        args.push_back(algo_opt);
    }
    args.push_back(dataname);
    if (repeat >= 0){
        args.push_back(to_string(repeat));
    }

    // BaseLine algorithm
    HSet *final = NULL;
    HSet *ret, *ret_nv;
    
    if(algorithm.compare("ns") == 0){
        cout << "Run NS" << endl;
        if (alpha >= 0.0){
            string alpha_str;
            stringstream stream;
            stream << std::fixed << std::setprecision(4) << alpha;
            alpha_str = stream.str();
            args[1] = args[1] + "_" + alpha_str;  
        }  
        outputdir = make_directory(args);
        cout << outputdir << endl;
        if (already_run(outputdir + "sampled_graph.txt", graph->number_of_hedges) && (!recalculate)){
            cout << "Already Run" << endl;
            final = NULL;
        }
        else{
            AlgorithmNS *algo = new AlgorithmNS(alpha, outputdir, algo_opt, graph);
            final = algo->run(accuracy);
            free(algo);
            algo = NULL;
        }
    }
    else if(((int)algorithm.size() >= 4) && (algorithm.compare("essz") == 0)){
        string alpha_str;
        stringstream stream;
        stream << std::fixed << std::setprecision(2) << alpha;
        alpha_str = stream.str();
        stream.str("");
        string beta_str;
        stream << std::fixed << std::setprecision(2) << beta;
        beta_str = stream.str();
        args[1] = args[1] + "_" + alpha_str + "_" + beta_str;
        outputdir = make_directory(args);
        cout << "Output directory is " << outputdir << endl;
        // run Sampling
        if (already_run(outputdir + "sampled_graph.txt", graph->number_of_hedges) && (!recalculate)){
            cout << "Already Run" << endl;
            final = NULL;
        }
        else{
            AlgorithmESSZ *algo = new AlgorithmESSZ(outputdir, algo_opt, alpha, beta, graph);
            final = algo->run(accuracy);
            free(algo);
            algo = NULL;
        }
    }
    else if(algorithm.compare("es") == 0){
        cout << "Run ES" << endl;
        if (alpha >= 0.0){
            string alpha_str;
            stringstream stream;
            stream << std::fixed << std::setprecision(4) << alpha;
            alpha_str = stream.str();
            args[1] = args[1] + "_" + alpha_str;
        }
        outputdir = make_directory(args);
        if (already_run(outputdir + "sampled_graph.txt", graph->number_of_hedges) && (!recalculate)){
            cout << "Already Run" << endl;
            final = NULL;
        }
        else{
            AlgorithmES *algo = new AlgorithmES(outputdir, algo_opt, alpha, graph);
            final = algo->run(accuracy);
            free(algo);
            algo = NULL;
        }
    }
    else if(algorithm.compare("ff") == 0){
        cout << "Run FF" << endl;
        stringstream stream;
        stream << std::fixed << std::setprecision(2) << p;
        string p_str = stream.str();

        stream.str("");
        stream << std::fixed << std::setprecision(2) << q;
        string q_str = stream.str();

        args[1] = args[1] + "_" + p_str + "_" + q_str;
        if (noinduce){
            args[1] = args[1] + "_noinduce";
        }
        outputdir = make_directory(args);
        if (already_run(outputdir + "sampled_graph.txt", graph->number_of_hedges) && (!recalculate)){
            cout << "Already Run" << endl;
            final = NULL;
        }
        else{
            Algorithm_FF *algo = new Algorithm_FF(p, q, algo_opt, outputdir, graph);
            final = algo->run(accuracy);
            free(algo);
            algo = NULL;
        }
    }
    else if(algorithm.compare("rw") == 0){
        cout << "Run RW" << endl;
        stringstream stream;
        stream << std::fixed << std::setprecision(2) << maxlength;
        string l_str = stream.str();
        stream.str("");
        stream << std::fixed << std::setprecision(2) << restart;
        string r_str = stream.str();

        args[1] = args[1] + "_" + l_str + "_" + r_str;
        if (noinduce){
            args[1] = args[1] + "_noinduce";
        }
        outputdir = make_directory(args);
        if (already_run(outputdir + "sampled_graph.txt", graph->number_of_hedges) && (!recalculate)){
            cout << "Already Run" << endl;
            final = NULL;
        }
        else{
            Algorithm_RW *algo = new Algorithm_RW(algo_opt, outputdir, graph, restart, maxlength, noinduce); 
            final = algo->run(accuracy);
            free(algo);
            algo = NULL;
        }
    }
    else if(algorithm.compare("tihs") == 0){
        cout << "Run TIHS" << endl;
        outputdir = make_directory(args);
        if (already_run(outputdir + "sampled_graph.txt", graph->number_of_hedges) && (!recalculate)){
            cout << "Already Run" << endl;
            final = NULL;
        }
        else{
            Algorithm_TIHS *algo = new Algorithm_TIHS(outputdir, graph); 
            final = algo->run(accuracy);
            free(algo);
            algo = NULL;
        }
    }

    // Other options
    if(algorithm.compare(0, 10, "helperdist") == 0){
        args.clear();
        args.push_back(inputdir);
        args.push_back(dataname);
        if (repeat >= 0){
            args.push_back(to_string(repeat));
        }
        outputdir = make_directory(args);
        cout << outputdir << endl;
        
        stringstream stream;
        string sp_str;
        stream << std::fixed << std::setprecision(1) << samplingportion;
        sp_str = stream.str();

        if ((int)algo_opt.size() == 0){
            algo_opt = "degree" + sp_str + ",size" + sp_str + ",pairdeg" + sp_str + ",intersection" + sp_str + ",wcc" + sp_str;
        }
        // prepare input
        vector<int> hyperedges;
        int unit_numhedge = (int)ceil((double)graph->number_of_hedges * (double)samplingportion);
        bool answer_flag = false;
        vector<string> tmp = split(outputdir, '/');
        for (int i = 0 ; i < (int)tmp.size() ; i++){
            if (tmp[i].compare("answer") == 0){
                answer_flag = true;
            }
        }
        if (answer_flag){
            for (int h = 0 ; h < graph->number_of_hedges; h++){
                hyperedges.push_back(h);
            }
        }
        else{
            if (!file_exist (outputdir +  "sampled_hindexes.txt")){
                return 0;
            }
            else if (!already_run(outputdir + "sampled_graph.txt", graph->number_of_hedges)){
                return 0;
            }
            string inputname = outputdir +  "sampled_hindexes.txt";
            ifstream inputFile(inputname.c_str());
            string line;
            while (getline(inputFile, line)){
                int h = stoi(line);
                hyperedges.push_back(h);
            }
        }
        if ((int)algo_opt.size() > 0){
            if (answer_flag){
                cout << "Calculate " <<  algo_opt << endl;
            }
            // remove outputs
            string output;
            tmp = split(algo_opt, ',');
            for (int i = 0 ; i < (int)tmp.size() ; i++){
                output = outputdir + tmp[i] + ".txt";
                remove(output.c_str());
            }
            // Evaluation
            int num_hedges = 0;
            set<int> init_set;
            HelperDist *algo = new HelperDist(init_set, graph, outputdir, algo_opt);
            cout << "Start Add" << endl;
            clock_t start = clock();
            for (int hi = 0; hi < graph->number_of_hedges ; hi++){
                int h = hyperedges[hi];
                init_set.insert(h);
                num_hedges++;
                if (num_hedges == unit_numhedge){
                    algo->update(init_set, graph);
                    algo->save_properties();
                    init_set.clear();
                    break;
                }
            }
        }
    }
    else if(algorithm.compare(0, 6, "helper") == 0){
        args.clear();
        args.push_back(inputdir);
        args.push_back(dataname);
        if (repeat >= 0){
            args.push_back(to_string(repeat));
        }
        outputdir = make_directory(args);
        cout << outputdir << endl;
        if ((int)algo_opt.size() == 0){
            algo_opt = "intersection,intersect_avg,densification,clusteringcoef,sizewcc";
        }
        // prepare input
        vector<int> hyperedges;
        int unit_numhedge = (int)ceil((double)graph->number_of_hedges / (double)accuracy);
        int num_step = 0;
        int num_hedges = 0;
        bool answer_flag = false;
        vector<string> tmp = split(outputdir, '/');
        for (int i = 0 ; i < (int)tmp.size() ; i++){
            if (tmp[i].compare("answer") == 0){
                answer_flag = true;
                for (int h = 0 ; h < graph->number_of_hedges; h++){
                    hyperedges.push_back(h);
                    num_hedges += 1;
                    if (( (num_hedges > 0) && (num_hedges % unit_numhedge == 0)) || (num_hedges == graph->number_of_hedges)){
                        num_step += 1;
                    }
                }
                break;
            }
        }
        if (!answer_flag){
            if (!file_exist (outputdir +  "sampled_hindexes.txt")){
                return 0;
            }
            else if (!already_run(outputdir + "sampled_graph.txt", graph->number_of_hedges)){
                return 0;
            }
            string inputname = outputdir +  "sampled_hindexes.txt";
            ifstream inputFile(inputname.c_str());
            string line;
            while (getline(inputFile, line)){
                int h = stoi(line);
                hyperedges.push_back(h);
                num_hedges += 1;
                if (( (num_hedges > 0) && (num_hedges % unit_numhedge == 0)) || (num_hedges == graph->number_of_hedges)){
                    num_step += 1;
                }
            }
        }
        // check already exist
        if (!recalculate){
            tmp = split(algo_opt, ',');
            string new_evalopt = "";
            for (int i = 0 ; i < (int)tmp.size() ; i++){
                if (!already_run(outputdir + tmp[i] + ".txt", num_step)){
                    if ((int)new_evalopt.size() == 0){
                        new_evalopt += tmp[i];
                    }
                    else{
                        new_evalopt += ',' + tmp[i];
                    }
                }
                else{
                    cout << "Already Run " << tmp[i] << endl;
                }
            }
            algo_opt = new_evalopt;
        }
        if ((int)algo_opt.size() > 0){
            if (answer_flag){
                cout << "Calculate " <<  algo_opt << endl;
            }
            // remove outputs
            string output;
            tmp = split(algo_opt, ',');
            for (int i = 0 ; i < (int)tmp.size() ; i++){
                output = outputdir + tmp[i] + ".txt";
                remove(output.c_str());
            }
            // Evaluation
            num_hedges = 0;
            set<int> init_set;
            Helper *algo = new Helper(init_set, graph, outputdir, algo_opt);
            cout << "Start Add" << endl;
            clock_t start = clock();
            for (int hi = 0; hi < graph->number_of_hedges ; hi++){
                int h = hyperedges[hi];
                init_set.insert(h);
                num_hedges++;
                if (( (num_hedges > 0) && (num_hedges % unit_numhedge == 0)) || (num_hedges == graph->number_of_hedges)){
                    int idx = (int)floor(num_hedges / unit_numhedge);
                    if ((idx % 10) == 0){
                        cout << "[" << idx << "]" << endl;
                        // cout << to_string((int)algo->neighbor.size()) << endl;
                        if (answer_flag){
                            clock_t end = clock();
                            double runtime =  ((double)(end - start) / CLOCKS_PER_SEC) / 60;
                            cout << "Time(min) = " << to_string(runtime) << endl;
                        }
                    }
                    algo->update(init_set, graph);
                    algo->save_properties();
                    init_set.clear();
                }
            }
        }
    }
    cout << "End Sampling" << endl;
    return 0;
}