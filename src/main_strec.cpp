/*
 * main_strec.cpp
 *
 */

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "state_reconstructor.h"
#include "rate_model.h"
#include "optimize_state_reconstructor_nlopt.h"
#include "optimize_state_reconstructor_periods_nlopt.h"
#include "tree_utils.h"
#include "log.h"


void print_help() {
    cout << "This will conduct state reconstruction analyses." << endl;
    cout << endl;
    cout << "Usage: pxstrec [OPTION]... [FILE]..." << endl;
    cout << endl;
    cout << " -d, --dataf=FILE    input data file" <<endl;
    cout << " -w, --datawide      data is in wide format so (001 instead of 2)" <<endl;
    cout << " -z, --dataz         data is in probability format (0,1,0,0)" << endl;
    cout << " -t, --treef=FILE    input tree file" << endl;
    cout << " -c, --conf=FILE     configuration file" << endl;
    cout << " -o, --outanc=FILE   output file for ancestral calc" << endl;
    cout << " -n, --outstnum=FILE output file for stochastic mapping number" << endl;
    cout << " -a, --outstnumany=FILE output file for stochastic mapping number any" << endl;
    cout << " -m, --outsttim=FILE output file for stochastic mapping duration" << endl;
    cout << " -p, --periods=NUMS  comma separated times" << endl;
    cout << " -l, --logf=FILE     log file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}
/*
 * add you name if you contribute (probably add another line)
 */
string versionline("pxstrec 0.2\nCopyright (C) 2017 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim), Joseph W. Brown");

static struct option const long_options[] =
{
    {"dataf", required_argument, NULL, 'd'},
    {"datawide", no_argument, NULL, 'w'},
    {"dataz", no_argument,  NULL, 'z'},
    {"treef", required_argument, NULL, 't'},
    {"conf", required_argument, NULL, 'c'},
    {"outanc", required_argument, NULL, 'o'},
    {"outstnum", required_argument, NULL, 'n'},
    {"outsttim", required_argument, NULL, 'm'},
    {"outstnumany", required_argument, NULL, 'a'},
    {"periods", required_argument, NULL, 'p'},
    {"logf", required_argument, NULL, 'l'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

bool checkdata(Tree * intree, vector<Sequence> runseqs);
bool checkdata(Tree * intree, vector<Sequence> runseqs) {
    vector<string> ret;
    set<string> seqnames;
    set<string> treenames;
    for (int i=0; i < intree->getExternalNodeCount(); i++) {
        treenames.insert(intree->getExternalNode(i)->getName());
    }
    for (unsigned int i=0; i < runseqs.size(); i++) {
        seqnames.insert(runseqs[i].get_id());
    }
    vector<string> v(treenames.size()+seqnames.size());
    vector<string>::iterator it;
    it=set_difference (seqnames.begin(), seqnames.end(), treenames.begin(), treenames.end(), v.begin());
    if (int(it-v.begin()) > 0) {
        cerr << "there are " << int(it - v.begin()) << " taxa that have the wrong names.\n";
    }
    for (it=v.begin() ; it != v.end(); it++) {
        if ((*it).size() > 1) {
            cout << *it << endl;
        }
    }
    return seqnames == treenames;
}

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool datafileset = false; // not used
    bool treefileset = false; // not used
    bool conffileset = false;
    bool logfileset = false;
    bool outancfileset = false;
    bool outstochtimefileset = false;
    bool outstochnumfileset = false;
    bool outstochnumanyfileset = false;
    bool datawide = false;
    bool periodsset = false;
    bool dataz = false; //the datafile will have probabilities
    char * conff = NULL;
    char * treef = NULL;
    char * dataf = NULL;
    char * logf = NULL;
    char * outanc = NULL;
    char * outnum = NULL;
    char * outnumany = NULL;
    char * outtime = NULL;
    string periodstring;
    vector<string> ptokens;
    vector<double> period_times;
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "d:t:c:o:n:m:a:l:p:hVwz", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'd':
                datafileset = true;
                dataf = strdup(optarg);
                check_file_exists(dataf);
                break;
            case 'z':
                dataz = true;
                datawide = true;
                break;
            case 'w':
                datawide = true;
                break;
            case 't':
                treefileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'c':
                conffileset = true;
                conff = strdup(optarg);
                check_file_exists(conff);
                break;
            case 'o':
                outancfileset = true;
                outanc = strdup(optarg);
                break;
            case 'n':
                outstochnumfileset = true;
                outnum = strdup(optarg);
                break;
            case 'm':
                outstochtimefileset = true;
                outtime = strdup(optarg);
                break;
            case 'a':
                outstochnumanyfileset = true;
                outnumany = strdup(optarg);
                break;
            case 'p':
                periodsset = true;
                periodstring = (strdup(optarg));
                parse_comma_list(periodstring, period_times);
                break;
            case 'l':
                logfileset = true;
                logf = strdup(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
            case 'V':
                cout << versionline << endl;
                exit(0);
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
    }
    
    ofstream * logout = NULL;
    ostream * loos = NULL;
    
    if (logfileset == true) {
        logout = new ofstream(logf);
        loos = logout;
    } else {
        loos = &cout;
    }

    if (conffileset == false) {
        cerr << "right now, you need to have a conf file. to change soon." << endl;
        print_help();
        exit(0);
    }

    bool verbose = true;
    string outfile_stochnum_any ="";
    string ratematrixfile = "";
    vector<vector<double> > ratematrix;
    bool estimate = true;
    map<string,vector<string> > mrcas;
    vector<string> stochtime;
    vector<string> stochnumber;
    vector<string> stochnumber_any;
    vector<string> ancstates;
    string freeparams = "_one_"; //right now, just _all_ or _one_

    /*************
     * read the configuration file
     **************/
    ifstream ifs(conff);
    string line;
    while (getline(ifs, line)) {
        if (line.size() > 1) {
            if (strcmp( (&line[0]), "#") != 0) {
            //if ((&line[0]) != "#") {
                vector<string> tokens;
                string del("=");
                tokens.clear();
                tokenize(line, tokens, del);
                for (unsigned int j=0; j < tokens.size(); j++) {
                    trim_spaces(tokens[j]);
                }
                if (!strcmp(tokens[0].c_str(), "freeparams")) {
                    freeparams = tokens[1];
                } else if (!strcmp(tokens[0].c_str(),  "outstnum_any")) {
                    outfile_stochnum_any = tokens[1];
                } else if (!strcmp(tokens[0].c_str(),  "ratematrix")) {
                    ratematrixfile = tokens[1];
                    if (ratematrixfile == "d" || ratematrixfile == "D") {
                        ratematrixfile = "";
                    }
                    estimate =false;
                } else if (!strcmp(tokens[0].c_str(), "mrca")) {
                    /*
                     * need to make sure that the mrca exist in the tree, 
                     * that the names are correct
                     */
                    vector<string> searchtokens;
                    tokenize(tokens[1], searchtokens, ",     ");
                    for (unsigned int j=0; j < searchtokens.size(); j++) {
                        trim_spaces(searchtokens[j]);
                    }
                    vector<string> mrc;
                    for (unsigned int j=1; j < searchtokens.size(); j++) {
                        mrc.push_back(searchtokens[j]);
                    }
                    mrcas[searchtokens[0]] = mrc;
                } else if (!strcmp(tokens[0].c_str(), "stochtime")) {
                    vector<string> searchtokens;
                    tokenize(tokens[1], searchtokens, ",     ");
                    for (unsigned int j=0; j < searchtokens.size(); j++) {
                        trim_spaces(searchtokens[j]);
                        stochtime.push_back(searchtokens[j]);
                    }
                } else if (!strcmp(tokens[0].c_str(), "stochnumber")) {
                    vector<string> searchtokens;
                    tokenize(tokens[1], searchtokens, ",     ");
                    for (unsigned int j=0; j < searchtokens.size(); j++) {
                        trim_spaces(searchtokens[j]);
                        stochnumber.push_back(searchtokens[j]);
                    }
                } else if (!strcmp(tokens[0].c_str(), "stochnumber_any")) {
                    vector<string> searchtokens;
                    tokenize(tokens[1], searchtokens, ",     ");
                    for (unsigned int j=0; j < searchtokens.size(); j++) {
                        trim_spaces(searchtokens[j]);
                        stochnumber_any.push_back(searchtokens[j]);
                    }
                } else if (!strcmp(tokens[0].c_str(), "ancstates")) {
                    vector<string> searchtokens;
                    tokenize(tokens[1], searchtokens, ",     ");
                    for (unsigned int j=0; j < searchtokens.size(); j++) {
                        trim_spaces(searchtokens[j]);
                        ancstates.push_back(searchtokens[j]);
                    }
                }
            }
        }
    }
    if (verbose) {
        (*loos) << "finished reading config file" << endl;
    }
    /**
     * read the data file
     */

    vector<Sequence> seqs;
    Sequence seq;
    ifstream * fstr = new ifstream(dataf);
    istream * pios = fstr;
    line = "";
    int ft = test_seq_filetype_stream(*pios,line);
    while (read_next_seq_from_stream(*pios,ft,line,seq)) {
    //(*loos) << seq.get_sequence() << endl;
        seqs.push_back(seq);
    }
    //fasta has a trailing one
    if (ft == 2) {
        seqs.push_back(seq);
    }
    if (verbose) {
        (*loos) << "taxa: " << seqs.size() << endl;
    }

    /**
     * read the tree file
     */
    TreeReader tr;
    vector<Tree *> trees;
    ifstream infile2(treef);
    if (!infile2) {
        cerr << "Could not open treefile." << endl;
        return 1;
    }
    line = "";
    while (getline(infile2, line)) {
        if (line.length() > 5) {
            trees.push_back(tr.readTree(line));
        }
    }
    infile2.close();
    if (verbose) {
        (*loos) << "trees: "<< trees.size() << endl;
    }
    
    /**
     * process the data
     * datawide means 
     * 
     * NOT datawide means
     */
    int nstates;
    int nsites;
    if (datawide) {
        if (dataz == false) {
            nstates = seqs[0].get_sequence().length();
        } else {
            vector<string> searchtokens;
            tokenize(seqs[0].get_sequence(), searchtokens, ",");
            nstates = searchtokens.size();
        }
        nsites = 1;
    } else {
        int maxstate=1;
        vector<string> searchtokens;
        tokenize(seqs[0].get_sequence(), searchtokens, " ");
        for (unsigned int j=0; j < searchtokens.size(); j++) {
            trim_spaces(searchtokens[j]);
        }
        nsites = searchtokens.size();
        if (verbose) {
            (*loos) << "nsites: " << nsites << endl;
        }
        for (unsigned int se = 0;se<seqs.size();se++) {
            searchtokens = vector<string> ();
            tokenize(seqs[se].get_sequence(), searchtokens, "     ");
            for (unsigned int j=0; j < searchtokens.size(); j++) {
                trim_spaces(searchtokens[j]);
                int pos = atoi(searchtokens[j].c_str());
                if (pos > maxstate) {
                    maxstate = pos;
                }
            }
        }
        nstates = maxstate+1;//TODO this can be determined by largest number +1
    }
    if (verbose) {
        (*loos) << "total number of states in dataset: " << nstates << endl;
    }
    //reading ratematrixfile
    if (ratematrixfile.size() > 1) {
        ratematrix = processRateMatrixConfigFile(ratematrixfile,nstates);
    }
    //end ratematrixfile

    ofstream ancout;
    ofstream stnumout;
    ofstream sttimeout;
    ofstream sttnumout_any;
    
    if (ancstates.size() > 0 && outancfileset == true) {
        ancout.open(outanc,ios::out);
        ancout << "site\ttree\tMRCA\tlnL";
        for (int i=0; i < nstates; i++) {
            ancout << "\tstate_" << i+1;
        }
        ancout << endl;
    }
    if (stochnumber.size() > 0 && outstochnumfileset == true) {
        stnumout.open(outnum,ios::out);
        stnumout << "site\ttree\tMRCA\tlnL";
        for (int i=0; i < nstates; i++) {
            for (int j=0; j < nstates; j++) {
                if (i != j) {
                    stnumout << "\tstate_" << i+1 << "->state_" << j+1;
                }
            }
        }
        stnumout << endl;
    }
    if (stochtime.size() > 0 && outstochtimefileset == true ) {
        sttimeout.open(outtime,ios::out);
        sttimeout << "site\ttree\tMRCA\tlnL";
        for (int i=0; i < nstates; i++) {
            sttimeout << "\tstate_" << i+1;
        }
        sttimeout << endl;
    }
    if (stochnumber_any.size() > 0 && outstochnumanyfileset == true) {
        sttnumout_any.open(outnumany,ios::out);
        sttnumout_any << "site\ttree\tMRCA\tlnL";
        sttnumout_any << "\tanystate";
        sttnumout_any << endl;
    }
    
    for (int n = 0; n < nsites; n++) {
        (*loos) << "site: " << n+1 << endl;
        /*
         * this converts the data and is a little long to accomodate datasets
         * with sites that don't have all the states but the results can still
         * be printed to the same outfile for analysis after
         */
        //need to put the data into a wide view
        vector<Sequence> runseqs;
        int nstates_site_n;
        vector<int> existing_states(nstates,0);
        if (datawide == false) {
            for (unsigned int se = 0;se<seqs.size();se++) {
                vector<string> searchtokens;
                tokenize(seqs[se].get_sequence(), searchtokens, "     ");
                for (unsigned int j=0; j < searchtokens.size(); j++) {
                    trim_spaces(searchtokens[j]);
                }
                string tseqs(nstates,'0');
                if (searchtokens[n]=="?") {
                    for (int mse = 0; mse < nstates; mse++) {
                    tseqs.replace(mse,1,"1");
                    }
                } else {
                    int pos = atoi(searchtokens[n].c_str());
                    tseqs.replace(pos,1,"1");
                }
                for (int i=0; i < nstates; i++) {
                    if (tseqs.at(i) =='1') {
                        existing_states[i] = 1;
                    }
                }
                Sequence tse =Sequence(seqs[se].get_id(),tseqs,true);
                runseqs.push_back(tse);
            }
            nstates_site_n = sum(existing_states);
        } else {
            runseqs = seqs;
            if (dataz == false) {
                for (unsigned int se=0; se < seqs.size(); se++) {
                    for (int i=0; i < nstates; i++) {
                        if (seqs[se].get_sequence().at(i) =='1') {
                            existing_states[i] = 1;
                        }
                    }
                }
                nstates_site_n = sum(existing_states);
            } else {
                for (int i=0; i < nstates; i++) {
                    existing_states[i] = 1;
                }
                nstates_site_n = nstates;
            }
        }
        //mapping the existing states to the full states
        //int statecnt = 0; // not used
        for (int i=nstates-1; i >= 0; i--) {
            if (existing_states[i] == 1) {
                continue;
            } else {
                for (unsigned int se=0; se < runseqs.size(); se++) {
                    runseqs[se].set_sequence(runseqs[se].get_sequence().erase(i,1));
                }
            }
        }

        if (verbose) {
            (*loos) <<"states: " << nstates_site_n << endl;
            (*loos) << "trees: ";
        }
        for (unsigned int i=0; i < trees.size(); i++) {
            if (verbose) {
                (*loos) << i << endl;
            }
            vector<RateModel> rms;
            RateModel rm(nstates_site_n);
            StateReconstructor sr(rm,rms);
            rm.setup_P(0.1,false);
            if (periodsset == true) {
                rms.push_back(rm);
                for (unsigned int p=1; p < period_times.size(); p++) {
                    RateModel rm2(nstates_site_n);
                    rm2.setup_P(0.1,false);
                    rms.push_back(rm2);
                }
                sr.set_periods(period_times,rms);
            }
            sr.set_store_p_matrices(false);
            Tree * tree = trees[i];
            if (verbose) {
                (*loos) << "tips: "<< tree->getExternalNodeCount() << endl;
            }
            sr.set_tree(tree);
            if (periodsset == true) {
                sr.set_periods_model();
            }
            //checking that the data and the tree have the same names
            if (checkdata(tree,runseqs) == 0) {
                exit(0);
            }
            bool same;
            if (dataz == false) {
                same = sr.set_tip_conditionals(runseqs);
            } else {
                same = sr.set_tip_conditionals_already_given(runseqs);
            }
            if (same == true) {
                (*loos) << "skipping calculation" <<endl;
                continue;
            }
            double finallike; Superdouble totlike_sd;
            if (periodsset == false) {
                mat free_var(nstates_site_n,nstates_site_n);
                free_var.fill(0);
                int ct = 0;
                if (freeparams == "_one_") {
                    ct = 1;
                } else if (freeparams == "_all_") {
                    ct = 0;
                    for (int k=0; k < nstates_site_n; k++) {
                        for (int j=0; j < nstates_site_n; j++) {
                            if (k != j) {
                                free_var(k,j) = ct;
                                ct += 1;
                            }
                        }
                    }
                }
                if (verbose) {
                    (*loos) << free_var << endl;
                    (*loos) << ct << endl;
                }
                rm.neg_p = false;
                cout << "likelihood: " << sr.eval_likelihood() << endl;
                //estimating the optimal rates
                if (estimate) {//optimize
                    optimize_sr_nlopt(&rm,&sr,&free_var,ct);
                } else { // requires that the ratematrix is available
                    for (int i=0; i < nstates_site_n; i++) {
                        for (int j=0; j < nstates_site_n; j++) {
                            free_var(i,j) = ratematrix[i][j];
                        }
                    }
                }
                //end estimating
                if (verbose) {
                    (*loos) << free_var << endl;
                }
                rm.setup_Q(free_var);
                sr.set_store_p_matrices(true);
                finallike = sr.eval_likelihood();
                if (verbose) {
                    (*loos) << "final_likelihood: " << finallike << endl;
                }
            } else { //optimize with periods
                vector<mat> periods_free_var(period_times.size());
                int ct = 0;
                if (freeparams == "_one_") {
                    ct = 1;
                    for (unsigned int s=0; s < period_times.size(); s++) {
                        mat free_var(nstates_site_n,nstates_site_n);free_var.fill(0);
                        periods_free_var[s] = free_var;
                    }
                } else if (freeparams == "_all_") {
                    ct = 0;
                    for (unsigned int s=0; s < period_times.size(); s++) {
                        mat free_var(nstates_site_n,nstates_site_n);free_var.fill(0);
                        for (int k=0; k < nstates_site_n; k++) {
                            for (int j=0; j < nstates_site_n; j++) {
                            if (k != j) {
                                free_var(k, j) = ct;
                                ct += 1;
                            }
                            }
                        }
                        periods_free_var[s] = free_var;
                    }
                }
                if (verbose) {
                    for (unsigned int s=0; s < period_times.size(); s++) {
                        (*loos) << periods_free_var[s] << endl;
                    }
                    (*loos) << ct << endl;
                }
                rm.neg_p = false;
                cout << "likelihood: " << sr.eval_likelihood() << endl;
                optimize_sr_periods_nlopt(&rms,&sr,&periods_free_var,ct);
                if (verbose) {
                    for (unsigned int s=0; s < period_times.size(); s++) {
                        (*loos) << periods_free_var[s] << endl;
                    }
                    (*loos) << ct << endl;
                    cout << "////////////////////////" << endl;
                }
                for (unsigned int s=0; s < period_times.size(); s++) {
                    rms[s].setup_Q(periods_free_var[s]);
                }
                sr.set_store_p_matrices(true);
                finallike = sr.eval_likelihood();
                if (verbose) {
                    (*loos) << "final_likelihood: " << finallike << endl;
                }
                cout << "period set and so no ancestral states just yet" << endl;
                continue;
            }
            if (verbose) {
                (*loos) << "ancestral states" <<endl;
            }
            sr.prepare_ancstate_reverse();
            for (unsigned int j=0; j < ancstates.size(); j++) {
            if (ancstates[j] == "_all_") {
                vector<Superdouble> lhoods;
                for (int l=0; l < tree->getInternalNodeCount(); l++) {
                    lhoods = sr.calculate_ancstate_reverse_sd(*tree->getInternalNode(l));
                    totlike_sd = calculate_vector_Superdouble_sum(lhoods);

                    //bool neg = false; // not used
                    int excount = 0;
                    double highest = 0;
                    int high = 0;
                    for (int k=0; k < nstates; k++) {
                        if (existing_states[k] == 1) {
                            if (double(lhoods[excount]/totlike_sd) > highest) {
                                highest= double(lhoods[excount]/totlike_sd);
                                high = k;
                            }
                            excount += 1;
                        }
                    }
                    std::string s;
                    std::stringstream out;
                    out << high;
                    tree->getInternalNode(l)->setName(out.str());
                }
                ancout << getNewickString(tree) << endl;
            } else {
                vector<Superdouble> lhoods;
                if (verbose) {
                    (*loos) <<"node: " << tree->getMRCA(mrcas[ancstates[j]])->getName() << "\tmrca: " << ancstates[j] <<  endl;
                }
                ancout << n+1 << "\t" << i+1 << "\t" << ancstates[j] << "\t" << finallike;
                lhoods = sr.calculate_ancstate_reverse_sd(*tree->getMRCA(mrcas[ancstates[j]]));
                totlike_sd = calculate_vector_Superdouble_sum(lhoods);
                bool neg = false;
                int excount = 0;
                for (int k=0; k < nstates; k++) {
                    if (existing_states[k] == 1) {
                        if (verbose) {
                            (*loos) << double(lhoods[excount]/totlike_sd) << " ";//"(" << lhoods[excount] << ") ";
                        }
                        ancout << "\t" << double(lhoods[excount]/totlike_sd);
                        if (double(lhoods[excount]/totlike_sd) < 0)
                        neg = true;
                        excount += 1;
                    } else {
                        if (verbose) {
                            (*loos) << "NA" << " ";
                            ancout << "\t" << "NA";
                        }
                    }
                }
                if (neg == true) {
                    exit(0);
                }
                ancout <<endl;
                if (verbose) {
                    (*loos) << endl;
                }
            }
            }
            if (verbose) {
                (*loos) << endl;
                (*loos) << "stochastic time" << endl;
            }

            for (unsigned int j=0; j < stochtime.size(); j++) {
            if (tree->getMRCA(mrcas[stochtime[j]])->isRoot() == false) {
                vector<double> lhoods;
                if (verbose) {
                    (*loos)  << "mrca: " << stochtime[j] <<  endl;
                }
                sttimeout << n+1 << "\t" << i+1 << "\t" << stochtime[j]<< "\t" << finallike;
                bool neg = false;
                int excount = 0;
                for (int k=0; k < nstates; k++) {
                    if (existing_states[k]==1) {
                        sr.prepare_stochmap_reverse_all_nodes(excount,excount);
                        sr.prepare_ancstate_reverse();
                        vector<double> stoch = sr.calculate_reverse_stochmap(*tree->getMRCA(mrcas[stochtime[j]]),true);
                        double tnum = sum(stoch)/double(totlike_sd);
                        double bl = tree->getMRCA(mrcas[stochtime[j]])->getBL();
                        if (verbose) {
                            (*loos) << tnum << " ";
                        }
                        sttimeout << "\t" << tnum/bl;
                        if (tnum < 0) {
                            neg = true;
                        }
                        excount += 1;
                    } else {
                        if (verbose) {
                            (*loos) << "NA" << " ";
                        }
                        sttimeout << "\t" << "NA";
                    }

                }
                sttimeout << endl;
                if (verbose) {
                    (*loos) << endl;
                }
                if (neg == true) {
                    exit(0);
                }
            }
            }
            if (verbose) {
                (*loos) << endl;
                (*loos) << "stochastic number" << endl;
            }
            for (unsigned int j=0; j < stochnumber.size(); j++) {
                if (tree->getMRCA(mrcas[stochnumber[j]])->isRoot() == false) {
                    vector<double> lhoods;
                    if (verbose) {
                        (*loos) << "mrca: " << stochnumber[j] <<  endl;
                    }
                    stnumout << n+1 << "\t" << i+1 << "\t" << stochnumber[j]<< "\t" << finallike;
                    bool neg = false;
                    int excount = 0;
                    for (int k=0; k < nstates; k++) {
                        if (existing_states[k]==1) {
                            int excount2 = 0;
                            for (int l=0; l < nstates; l++) {
                                if (existing_states[l] == 1) {
                                    if (k == l) {
                                        if (verbose) {
                                            (*loos) << " - ";
                                        }
                                    } else {
                                        sr.prepare_stochmap_reverse_all_nodes(excount,excount2);
                                        sr.prepare_ancstate_reverse();
                                        vector<double> stoch = sr.calculate_reverse_stochmap(*tree->getMRCA(mrcas[stochnumber[j]]),false);
                                        double tnum = sum(stoch)/totlike_sd;
                                        if (verbose) {
                                            (*loos) << tnum << " ";
                                        }
                                        stnumout << "\t" << tnum;
                                        if (tnum < 0) {
                                            neg = true;
                                        }
                                    }
                                    excount2 += 1;
                                } else {
                                    if (verbose) {
                                        (*loos) << "NA" << " ";
                                    }
                                    stnumout << "\t" << "NA";
                                }
                            }
                            if (verbose) {
                                (*loos) << endl;
                            }
                            excount += 1;
                        } else {
                            for (int l=0; l < nstates; l++) {
                                if (k == l) {
                                    if (verbose) {
                                        (*loos) << " - ";
                                    }
                                } else {
                                    if (verbose) {
                                        (*loos) << "NA" << " ";
                                    }
                                    stnumout << "\t" << "NA";
                                }
                            }
                            if (verbose) {
                                (*loos) << endl;
                            }
                        }
                    }
                    stnumout << endl;
                    if (verbose) {
                        (*loos) << endl;
                    }
                    if (neg == true) {
                        exit(0);
                    }
                }
            }
            if (verbose) {
                (*loos) << endl;
            }
            if (verbose) {
                (*loos) << "stochastic number (any)" << endl;
            }
            if (stochnumber_any.size() > 0) {
                sr.prepare_stochmap_reverse_all_nodes_all_matrices();
                sr.prepare_ancstate_reverse();
            }
            for (unsigned int j=0; j < stochnumber_any.size(); j++) {
                if (tree->getMRCA(mrcas[stochnumber_any[j]])->isRoot() == false) {
                    vector<double> lhoods;
                    if (verbose) {
                        (*loos) <<"node: " << tree->getMRCA(mrcas[stochnumber_any[j]])->getName() << " mrca: " << stochnumber_any[j] <<  endl;
                    }
                    sttnumout_any << n+1 << "\t" << i+1 << "\t" << stochnumber_any[j]<< "\t" << finallike;
                    vector<double> stoch = sr.calculate_reverse_stochmap(*tree->getMRCA(mrcas[stochnumber_any[j]]),false);
                    double tnum = sum(stoch)/totlike_sd;
                    //(*loos) << sum(stoch)<< " "<<totlike << endl;
                    if (verbose) {
                        (*loos) << tnum << " " ;
                    }
                    sttnumout_any << "\t" << tnum;
                    sttnumout_any << endl;
                    if (verbose) {
                        (*loos) << endl;
                    }
                }
            }

            //delete tree;
        }
    }
    if (ancstates.size() > 0  && outancfileset == true) {
        ancout.close();
    }
    if (stochnumber.size() > 0) {
        stnumout.close();
    }
    if (stochtime.size() > 0) {
        sttimeout.close();
    }
    
    if (logfileset) {
        logout->close();
        delete loos;
    }
    return EXIT_SUCCESS;
}
