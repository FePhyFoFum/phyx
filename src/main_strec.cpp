#include <getopt.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>

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
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << "This will conduct state reconstruction analyses." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxstrec [OPTIONS]... FILES" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -d, --dataf=FILE    input data file" << std::endl;
    std::cout << " -w, --datawide      data is in wide format so (001 instead of 2)" << std::endl;
    std::cout << " -z, --dataz         data is in probability format (0,1,0,0)" << std::endl;
    std::cout << " -t, --treef=FILE    input tree file" << std::endl;
    std::cout << " -c, --conf=FILE     configuration file" << std::endl;
    std::cout << " -n, --outstnum=FILE output file for stochastic mapping number" << std::endl;
    std::cout << " -a, --outstnumany=FILE output file for stochastic mapping number any" << std::endl;
    std::cout << " -m, --outsttim=FILE output file for stochastic mapping duration" << std::endl;
    std::cout << " -p, --periods=NUMS  comma separated times" << std::endl;
    std::cout << " -l, --logf=FILE     log file, STOUT otherwise" << std::endl;
    std::cout << " -o, --outanc=FILE   output file for ancestral calc" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxstrec 1.2\n";
    vl += "Copyright (C) 2017-2021 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"dataf", required_argument, nullptr, 'd'},
    {"datawide", no_argument, nullptr, 'w'},
    {"dataz", no_argument,  nullptr, 'z'},
    {"treef", required_argument, nullptr, 't'},
    {"conf", required_argument, nullptr, 'c'},
    {"outanc", required_argument, nullptr, 'o'},
    {"outstnum", required_argument, nullptr, 'n'},
    {"outsttim", required_argument, nullptr, 'm'},
    {"outstnumany", required_argument, nullptr, 'a'},
    {"periods", required_argument, nullptr, 'p'},
    {"logf", required_argument, nullptr, 'l'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

bool checkdata(Tree * intree, std::vector<Sequence> runseqs);
bool checkdata(Tree * intree, std::vector<Sequence> runseqs) {
    std::vector<std::string> ret;
    std::set<std::string> seqnames;
    std::set<std::string> treenames;
    for (int i = 0; i < intree->getExternalNodeCount(); i++) {
        treenames.insert(intree->getExternalNode(i)->getName());
    }
    for (unsigned int i = 0; i < runseqs.size(); i++) {
        seqnames.insert(runseqs[i].get_id());
    }
    std::vector<std::string> v(treenames.size()+seqnames.size());
    std::vector<std::string>::iterator it;
    it = set_difference (seqnames.begin(), seqnames.end(), treenames.begin(), treenames.end(), v.begin());
    if (int(it-v.begin()) > 0) {
        std::cerr << "Error: there are " << int(it - v.begin())
                << " taxa that have the wrong names. Exiting." << std::endl;
        exit(1);
    }
    for (it = v.begin(); it != v.end(); it++) {
        if ((*it).size() > 1) {
            std::cout << *it << std::endl;
        }
    }
    return seqnames == treenames;
}

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    //bool datafileset = false; // not used
    //bool treefileset = false; // not used
    bool conffileset = false;
    bool logfileset = false;
    bool outancfileset = false;
    bool outstochtimefileset = false;
    bool outstochnumfileset = false;
    bool outstochnumanyfileset = false;
    bool datawide = false;
    bool periodsset = false;
    bool dataz = false; //the datafile will have probabilities
    char * conff = nullptr;
    char * treef = nullptr;
    char * dataf = nullptr;
    char * logf = nullptr;
    char * outanc = nullptr;
    char * outnum = nullptr;
    char * outnumany = nullptr;
    char * outtime = nullptr;
    std::string periodstring;
    std::vector<std::string> ptokens;
    std::vector<double> period_times;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "d:t:c:o:n:m:a:l:p:hVwzC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'd':
                //datafileset = true;
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
                //treefileset = true;
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
                std::cout << get_version_line() << std::endl;
                exit(0);
            case 'C':
                std::cout << get_phyx_citation() << std::endl;
                exit(0);
            default:
                print_error(argv[0]);
                exit(0);
        }
    }
    
    std::ofstream * logout = nullptr;
    std::ostream * loos = nullptr;
    
    if (logfileset) {
        logout = new std::ofstream(logf);
        loos = logout;
    } else {
        loos = &cout;
    }

    if (!conffileset) {
        std::cerr << "Error: you need to have a conf file. May change soon. Exiting." << std::endl;
        print_help();
        exit(0);
    }

    bool verbose = true;
    std::string outfile_stochnum_any;
    std::string ratematrixfile;
    std::vector<std::vector<double> > ratematrix;
    bool estimate = true;
    std::map<std::string, std::vector<std::string> > mrcas;
    std::vector<std::string> stochtime;
    std::vector<std::string> stochnumber;
    std::vector<std::string> stochnumber_any;
    std::vector<std::string> ancstates;
    std::string freeparams = "_one_"; //right now, just _all_ or _one_

    /*************
     * read the configuration file
     **************/
    std::ifstream ifs(conff);
    std::string line;
    std::vector<std::string> tokens;
    while (getline_safe(ifs, line)) {
        if (line.size() > 1) {
            if (strcmp( (&line[0]), "#") != 0) {
            //if ((&line[0]) != "#") {
                std::string del("=");
                tokens.clear();
                tokenize(line, tokens, del);
                for (unsigned int j = 0; j < tokens.size(); j++) {
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
                    std::vector<std::string> searchtokens;
                    tokenize(tokens[1], searchtokens, ",     ");
                    for (unsigned int j = 0; j < searchtokens.size(); j++) {
                        trim_spaces(searchtokens[j]);
                    }
                    std::vector<std::string> mrc;
                    for (unsigned int j = 1; j < searchtokens.size(); j++) {
                        mrc.push_back(searchtokens[j]);
                    }
                    mrcas[searchtokens[0]] = mrc;
                } else if (!strcmp(tokens[0].c_str(), "stochtime")) {
                    std::vector<std::string> searchtokens;
                    tokenize(tokens[1], searchtokens, ",     ");
                    for (unsigned int j = 0; j < searchtokens.size(); j++) {
                        trim_spaces(searchtokens[j]);
                        stochtime.push_back(searchtokens[j]);
                    }
                } else if (!strcmp(tokens[0].c_str(), "stochnumber")) {
                    std::vector<std::string> searchtokens;
                    tokenize(tokens[1], searchtokens, ",     ");
                    for (unsigned int j = 0; j < searchtokens.size(); j++) {
                        trim_spaces(searchtokens[j]);
                        stochnumber.push_back(searchtokens[j]);
                    }
                } else if (!strcmp(tokens[0].c_str(), "stochnumber_any")) {
                    std::vector<std::string> searchtokens;
                    tokenize(tokens[1], searchtokens, ",     ");
                    for (unsigned int j = 0; j < searchtokens.size(); j++) {
                        trim_spaces(searchtokens[j]);
                        stochnumber_any.push_back(searchtokens[j]);
                    }
                } else if (!strcmp(tokens[0].c_str(), "ancstates")) {
                    std::vector<std::string> searchtokens;
                    tokenize(tokens[1], searchtokens, ",     ");
                    for (unsigned int j = 0; j < searchtokens.size(); j++) {
                        trim_spaces(searchtokens[j]);
                        ancstates.push_back(searchtokens[j]);
                    }
                }
            }
        }
    }
    if (verbose) {
        (*loos) << "finished reading config file" << std::endl;
    }
    /**
     * read the data file
     */

    std::vector<Sequence> seqs;
    Sequence seq;
    std::ifstream * fstr = new std::ifstream(dataf);
    std::istream * pios = fstr;
    line = "";
    int ft = test_seq_filetype_stream(*pios, line);
    while (read_next_seq_from_stream(*pios, ft, line, seq)) {
    //(*loos) << seq.get_sequence() << std::endl;
        seqs.push_back(seq);
    }
    //fasta has a trailing one
    if (ft == 2) {
        seqs.push_back(seq);
    }
    if (verbose) {
        (*loos) << "taxa: " << seqs.size() << std::endl;
    }

    /**
     * read the tree file
     */
    TreeReader tr;
    std::vector<Tree *> trees;
    std::ifstream infile2(treef);
    if (!infile2) {
        std::cerr << "Error: could not open treefile. Exiting." << std::endl;
        exit(1);
    }
    line = "";
    while (getline_safe(infile2, line)) {
        if (line.length() > 5) {
            trees.push_back(tr.readTree(line));
        }
    }
    infile2.close();
    if (verbose) {
        (*loos) << "trees: " << trees.size() << std::endl;
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
        if (!dataz) {
            nstates = seqs[0].get_length();
        } else {
            std::vector<std::string> searchtokens;
            tokenize(seqs[0].get_sequence(), searchtokens, ",");
            nstates = searchtokens.size();
        }
        nsites = 1;
    } else {
        int maxstate=1;
        std::vector<std::string> searchtokens;
        tokenize(seqs[0].get_sequence(), searchtokens, " ");
        for (unsigned int j = 0; j < searchtokens.size(); j++) {
            trim_spaces(searchtokens[j]);
        }
        nsites = searchtokens.size();
        if (verbose) {
            (*loos) << "nsites: " << nsites << std::endl;
        }
        for (unsigned int se = 0;se<seqs.size();se++) {
            searchtokens = std::vector<std::string> ();
            tokenize(seqs[se].get_sequence(), searchtokens, "     ");
            for (unsigned int j = 0; j < searchtokens.size(); j++) {
                trim_spaces(searchtokens[j]);
                int pos = std::atoi(searchtokens[j].c_str());
                if (pos > maxstate) {
                    maxstate = pos;
                }
            }
        }
        nstates = maxstate+1;//TODO this can be determined by largest number +1
    }
    if (verbose) {
        (*loos) << "total number of states in dataset: " << nstates << std::endl;
    }
    //reading ratematrixfile
    if (ratematrixfile.size() > 1) {
        ratematrix = processRateMatrixConfigFile(ratematrixfile, nstates);
    }
    //end ratematrixfile
    std::ofstream ancout;
    std::ofstream stnumout;
    std::ofstream sttimeout;
    std::ofstream sttnumout_any;
    
    if (!ancstates.empty() && outancfileset) {
        ancout.open(outanc, std::ios::out);
        ancout << "site\ttree\tMRCA\tlnL";
        for (int i = 0; i < nstates; i++) {
            ancout << "\tstate_" << i+1;
        }
        ancout << std::endl;
    }
    if (!stochnumber.empty() && outstochnumfileset) {
        stnumout.open(outnum, std::ios::out);
        stnumout << "site\ttree\tMRCA\tlnL";
        for (int i = 0; i < nstates; i++) {
            for (int j = 0; j < nstates; j++) {
                if (i != j) {
                    stnumout << "\tstate_" << i+1 << "->state_" << j+1;
                }
            }
        }
        stnumout << std::endl;
    }
    if (!stochtime.empty() && outstochtimefileset ) {
        sttimeout.open(outtime, std::ios::out);
        sttimeout << "site\ttree\tMRCA\tlnL";
        for (int i = 0; i < nstates; i++) {
            sttimeout << "\tstate_" << i+1;
        }
        sttimeout << std::endl;
    }
    if (!stochnumber_any.empty() && outstochnumanyfileset) {
        sttnumout_any.open(outnumany, std::ios::out);
        sttnumout_any << "site\ttree\tMRCA\tlnL";
        sttnumout_any << "\tanystate";
        sttnumout_any << std::endl;
    }
    
    for (int n = 0; n < nsites; n++) {
        (*loos) << "site: " << n+1 << std::endl;
        /*
         * this converts the data and is a little long to accomodate datasets
         * with sites that don't have all the states but the results can still
         * be printed to the same outfile for analysis after
         */
        //need to put the data into a wide view
        std::vector<Sequence> runseqs;
        int nstates_site_n;
        std::vector<int> existing_states(nstates, 0);
        if (!datawide) {
            for (unsigned int se = 0;se<seqs.size();se++) {
                std::vector<std::string> searchtokens;
                tokenize(seqs[se].get_sequence(), searchtokens, "     ");
                for (unsigned int j = 0; j < searchtokens.size(); j++) {
                    trim_spaces(searchtokens[j]);
                }
                std::string tseqs(nstates, '0');
                if (searchtokens[n]=="?") {
                    for (int mse = 0; mse < nstates; mse++) {
                        tseqs.replace(mse, 1, "1");
                    }
                } else {
                    int pos = std::atoi(searchtokens[n].c_str());
                    tseqs.replace(pos, 1, "1");
                }
                for (int i = 0; i < nstates; i++) {
                    if (tseqs.at(i) == '1') {
                        existing_states[i] = 1;
                    }
                }
                Sequence tse = Sequence(seqs[se].get_id(), tseqs, true);
                runseqs.push_back(tse);
            }
            nstates_site_n = sum(existing_states);
        } else {
            runseqs = seqs;
            if (!dataz) {
                for (unsigned int se=0; se < seqs.size(); se++) {
                    for (int i = 0; i < nstates; i++) {
                        if (seqs[se].get_sequence().at(i) == '1') {
                            existing_states[i] = 1;
                        }
                    }
                }
                nstates_site_n = sum(existing_states);
            } else {
                for (int i = 0; i < nstates; i++) {
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
                    runseqs[se].set_sequence(runseqs[se].get_sequence().erase(i, 1));
                }
            }
        }

        if (verbose) {
            (*loos) << "states: " << nstates_site_n << std::endl;
            (*loos) << "trees: ";
        }
        for (unsigned int i = 0; i < trees.size(); i++) {
            if (verbose) {
                (*loos) << i << std::endl;
            }
            std::vector<RateModel> rms;
            RateModel rm(nstates_site_n);
            StateReconstructor sr(rm, rms);
            rm.setup_P(0.1, false);
            if (periodsset) {
                rms.push_back(rm);
                for (unsigned int p=1; p < period_times.size(); p++) {
                    RateModel rm2(nstates_site_n);
                    rm2.setup_P(0.1, false);
                    rms.push_back(rm2);
                }
                sr.set_periods(period_times, rms);
            }
            sr.set_store_p_matrices(false);
            Tree * tree = trees[i];
            if (verbose) {
                (*loos) << "tips: " << tree->getExternalNodeCount() << std::endl;
            }
            sr.set_tree(tree);
            if (periodsset) {
                sr.set_periods_model();
            }
            //checking that the data and the tree have the same names
            if (checkdata(tree, runseqs) == 0) {
                exit(0);
            }
            bool same;
            if (!dataz) {
                same = sr.set_tip_conditionals(runseqs);
            } else {
                same = sr.set_tip_conditionals_already_given(runseqs);
            }
            if (same) {
                (*loos) << "skipping calculation" << std::endl;
                continue;
            }
            double finallike; Superdouble totlike_sd;
            if (!periodsset) {
                mat free_var(nstates_site_n, nstates_site_n);
                free_var.fill(0);
                int ct = 0;
                if (freeparams == "_one_") {
                    ct = 1;
                } else if (freeparams == "_all_") {
                    ct = 0;
                    for (int k = 0; k < nstates_site_n; k++) {
                        for (int j = 0; j < nstates_site_n; j++) {
                            if (k != j) {
                                free_var(k, j) = ct;
                                ct += 1;
                            }
                        }
                    }
                }
                if (verbose) {
                    (*loos) << free_var << std::endl;
                    (*loos) << ct << std::endl;
                }
                rm.neg_p = false;
                std::cout << "likelihood: " << sr.eval_likelihood() << std::endl;
                //estimating the optimal rates
                if (estimate) {//optimize
                    optimize_sr_nlopt(&rm, &sr, &free_var, ct);
                } else { // requires that the ratematrix is available
                    for (int i = 0; i < nstates_site_n; i++) {
                        for (int j = 0; j < nstates_site_n; j++) {
                            free_var(i, j) = ratematrix[i][j];
                        }
                    }
                }
                //end estimating
                if (verbose) {
                    (*loos) << free_var << std::endl;
                }
                rm.setup_Q(free_var);
                sr.set_store_p_matrices(true);
                finallike = sr.eval_likelihood();
                if (verbose) {
                    (*loos) << "final_likelihood: " << finallike << std::endl;
                }
            } else { //optimize with periods
                std::vector<mat> periods_free_var(period_times.size());
                int ct = 0;
                if (freeparams == "_one_") {
                    ct = 1;
                    for (unsigned int s=0; s < period_times.size(); s++) {
                        mat free_var(nstates_site_n, nstates_site_n);
                        free_var.fill(0);
                        periods_free_var[s] = free_var;
                    }
                } else if (freeparams == "_all_") {
                    ct = 0;
                    for (unsigned int s=0; s < period_times.size(); s++) {
                        mat free_var(nstates_site_n, nstates_site_n);
                        free_var.fill(0);
                        for (int k = 0; k < nstates_site_n; k++) {
                            for (int j = 0; j < nstates_site_n; j++) {
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
                        (*loos) << periods_free_var[s] << std::endl;
                    }
                    (*loos) << ct << std::endl;
                }
                rm.neg_p = false;
                std::cout << "likelihood: " << sr.eval_likelihood() << std::endl;
                optimize_sr_periods_nlopt(&rms, &sr, &periods_free_var, ct);
                if (verbose) {
                    for (unsigned int s=0; s < period_times.size(); s++) {
                        (*loos) << periods_free_var[s] << std::endl;
                    }
                    (*loos) << ct << std::endl;
                    std::cout << "////////////////////////" << std::endl;
                }
                for (unsigned int s=0; s < period_times.size(); s++) {
                    rms[s].setup_Q(periods_free_var[s]);
                }
                sr.set_store_p_matrices(true);
                finallike = sr.eval_likelihood();
                if (verbose) {
                    (*loos) << "final_likelihood: " << finallike << std::endl;
                }
                std::cout << "period set and so no ancestral states just yet" << std::endl;
                continue;
            }
            if (verbose) {
                (*loos) << "ancestral states" << std::endl;
            }
            sr.prepare_ancstate_reverse();
            for (unsigned int j = 0; j < ancstates.size(); j++) {
            if (ancstates[j] == "_all_") {
                std::vector<Superdouble> lhoods;
                for (int l = 0; l < tree->getInternalNodeCount(); l++) {
                    lhoods = sr.calculate_ancstate_reverse_sd(*tree->getInternalNode(l));
                    totlike_sd = calculate_vector_Superdouble_sum(lhoods);

                    //bool neg = false; // not used
                    int excount = 0;
                    double highest = 0;
                    int high = 0;
                    for (int k = 0; k < nstates; k++) {
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
                ancout << getNewickString(tree) << std::endl;
            } else {
                std::vector<Superdouble> lhoods;
                if (verbose) {
                    (*loos) << "node: " << tree->getMRCA(mrcas[ancstates[j]])->getName()
                            << "\tmrca: " << ancstates[j] << std::endl;
                }
                ancout << n+1 << "\t" << i+1 << "\t" << ancstates[j] << "\t" << finallike;
                lhoods = sr.calculate_ancstate_reverse_sd(*tree->getMRCA(mrcas[ancstates[j]]));
                totlike_sd = calculate_vector_Superdouble_sum(lhoods);
                bool neg = false;
                int excount = 0;
                for (int k = 0; k < nstates; k++) {
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
                if (neg) {
                    exit(0);
                }
                ancout << std::endl;
                if (verbose) {
                    (*loos) << std::endl;
                }
            }
            }
            if (verbose) {
                (*loos) << std::endl;
                (*loos) << "stochastic time" << std::endl;
            }

            for (unsigned int j = 0; j < stochtime.size(); j++) {
            if (!tree->getMRCA(mrcas[stochtime[j]])->isRoot()) {
                std::vector<double> lhoods;
                if (verbose) {
                    (*loos)  << "mrca: " << stochtime[j] << std::endl;
                }
                sttimeout << n+1 << "\t" << i+1 << "\t" << stochtime[j]<< "\t" << finallike;
                bool neg = false;
                int excount = 0;
                for (int k = 0; k < nstates; k++) {
                    if (existing_states[k]==1) {
                        sr.prepare_stochmap_reverse_all_nodes(excount, excount);
                        sr.prepare_ancstate_reverse();
                        std::vector<double> stoch;
                        stoch = sr.calculate_reverse_stochmap(*tree->getMRCA(mrcas[stochtime[j]]), true);
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
                sttimeout << std::endl;
                if (verbose) {
                    (*loos) << std::endl;
                }
                if (neg) {
                    exit(0);
                }
            }
            }
            if (verbose) {
                (*loos) << std::endl;
                (*loos) << "stochastic number" << std::endl;
            }
            for (unsigned int j = 0; j < stochnumber.size(); j++) {
                if (!tree->getMRCA(mrcas[stochnumber[j]])->isRoot()) {
                    std::vector<double> lhoods;
                    if (verbose) {
                        (*loos) << "mrca: " << stochnumber[j] << std::endl;
                    }
                    stnumout << n+1 << "\t" << i+1 << "\t" << stochnumber[j]<< "\t" << finallike;
                    bool neg = false;
                    int excount = 0;
                    for (int k = 0; k < nstates; k++) {
                        if (existing_states[k]==1) {
                            int excount2 = 0;
                            for (int l = 0; l < nstates; l++) {
                                if (existing_states[l] == 1) {
                                    if (k == l) {
                                        if (verbose) {
                                            (*loos) << " - ";
                                        }
                                    } else {
                                        sr.prepare_stochmap_reverse_all_nodes(excount, excount2);
                                        sr.prepare_ancstate_reverse();
                                        std::vector<double> stoch;
                                        stoch = sr.calculate_reverse_stochmap(*tree->getMRCA(mrcas[stochnumber[j]]), false);
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
                                (*loos) << std::endl;
                            }
                            excount += 1;
                        } else {
                            for (int l = 0; l < nstates; l++) {
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
                                (*loos) << std::endl;
                            }
                        }
                    }
                    stnumout << std::endl;
                    if (verbose) {
                        (*loos) << std::endl;
                    }
                    if (neg) {
                        exit(0);
                    }
                }
            }
            if (verbose) {
                (*loos) << std::endl;
            }
            if (verbose) {
                (*loos) << "stochastic number (any)" << std::endl;
            }
            if (!stochnumber_any.empty()) {
                sr.prepare_stochmap_reverse_all_nodes_all_matrices();
                sr.prepare_ancstate_reverse();
            }
            for (unsigned int j = 0; j < stochnumber_any.size(); j++) {
                if (!tree->getMRCA(mrcas[stochnumber_any[j]])->isRoot()) {
                    //std::vector<double> lhoods; // not used
                    std::vector<double> stoch;
                    if (verbose) {
                        (*loos) << "node: " << tree->getMRCA(mrcas[stochnumber_any[j]])->getName()
                                << " mrca: " << stochnumber_any[j] << std::endl;
                    }
                    sttnumout_any << n+1 << "\t" << i+1 << "\t" << stochnumber_any[j]<< "\t" << finallike;
                    stoch = sr.calculate_reverse_stochmap(*tree->getMRCA(mrcas[stochnumber_any[j]]), false);
                    double tnum = sum(stoch)/totlike_sd;
                    //(*loos) << sum(stoch) << " " << totlike << std::endl;
                    if (verbose) {
                        (*loos) << tnum << " ";
                    }
                    sttnumout_any << "\t" << tnum;
                    sttnumout_any << std::endl;
                    if (verbose) {
                        (*loos) << std::endl;
                    }
                }
            }

            //delete tree;
        }
    }
    if (!ancstates.empty() && outancfileset) {
        ancout.close();
    }
    if (!stochnumber.empty()) {
        stnumout.close();
    }
    if (!stochtime.empty()) {
        sttimeout.close();
    }
    
    if (logfileset) {
        logout->close();
        delete loos;
    }
    return EXIT_SUCCESS;
}
