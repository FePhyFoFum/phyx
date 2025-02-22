#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <getopt.h>

#include "utils.h"
#include "log_manip.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "MCMC log file manipulator." << std::endl;
    std::cout << "Can combine and resample parameters or trees across files." << std::endl;
    std::cout << "Log files need not contain the same number of samples." << std::endl;
    std::cout << "Input files may be indicated using wildcards e.g. '*.trees'" << std::endl;
    std::cout << "Parameter log files are expected to be whitespace delimited." << std::endl;
    std::cout << "*NOTE* All values are in terms of number of SAMPLES (NOT generations)." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxlog [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -p, --parmf=FILE    input parameter log file(s)" << std::endl;
    std::cout << " -t, --treef=FILE    input tree log file(s)" << std::endl;
    std::cout << " -b, --burnin=INT    number of samples to exclude at the beginning of a file" << std::endl;
    std::cout << " -n, --thin=INT      interval of resampling" << std::endl;
    // not yet implemented:
//    std::cout << " -r, --rand=INT      number of random samples (without replacement)" << std::endl;
    std::cout << " -i, --info          calculate log file attributes and exit" << std::endl;
    std::cout << " -s, --summarize     summary statistics of samples (parameter logs only)" << std::endl;
    std::cout << " -c, --columns       print out column names (parameter logs only)" << std::endl;
    std::cout << " -d, --delete=CSL    delete parm columns by 1-index sep by commas (NO SPACES!)" << std::endl;
    std::cout << " -k, --keep=CSL      keep only parm columns by 1-index sep by commas (NO SPACES!)" << std::endl;
//    std::cout << " -x, --seed=INT      random number seed, clock otherwise" << std::endl;
    std::cout << " -v, --verbose       make the output more verbose" << std::endl;
    std::cout << " -o, --outf=FILE     output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxlog 1.3.2\n";
    vl += "Copyright (C) 2016-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"parmf", required_argument, nullptr, 'p'},
    {"treef", required_argument, nullptr, 't'},
    {"outf", required_argument, nullptr, 'o'},
    {"burnin", required_argument, nullptr, 'b'},
    {"thin", required_argument, nullptr, 'n'},
    {"rand", required_argument, nullptr, 'r'},
    {"info", no_argument, nullptr, 'i'},
    {"summarize", no_argument, nullptr, 's'},
    {"columns", no_argument, nullptr, 'c'},
    {"delete", required_argument, nullptr, 'd'},
    {"keep", required_argument, nullptr, 'k'},
    {"seed", required_argument, nullptr, 'x'},
    {"verbose", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool pfileset = false;
    bool tfileset = false;
    std::vector<std::string> input_files;
    std::vector<int> col_indices;
    
    int burnin = 0;
    int nthin = 1;
    int nrandom = -1;
    long int seed = -1;
    bool verbose = false;
    bool count = false;
    bool summarize = false;
    bool get_columns = false;
    bool delete_columns = false;
    bool keep_columns = false;
    std::string incolids;
    std::string logtype;
    char * outf = nullptr;
    bool argspresent = false;
    
    while (true) {
        int oi = -1;
        int curind;
        int c = getopt_long(argc, argv, "p:t:o:b:n:r:iscd:k:x:vhVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'p':
                pfileset = true;
                curind = optind - 1;
                while (curind < argc) {
                    std::string temp = strdup(argv[curind]);
                    curind++;
                    if (temp[0] != '-') {
                        std::ifstream infile(temp.c_str());
                        if (infile.good()) { // check that file exists
                            input_files.push_back(temp);
                            infile.close();
                        } else {
                            std::cerr << "Error: cannot find input file '" << temp << "'. Exiting." << std::endl;
                            exit(0);
                        }
                    } else {
                        optind = curind - 1;
                        break;
                    }
                }
                logtype = "parameter";
                argspresent = true;
                break;
            case 't':
                tfileset = true;
                curind = optind - 1;
                while (curind < argc) {
                    std::string temp = strdup(argv[curind]);
                    curind++;
                    if (temp[0] != '-') {
                        std::ifstream infile(temp.c_str());
                        if (infile.good()) { // check that file exists
                            input_files.push_back(temp);
                            infile.close();
                        } else {
                            std::cerr << "Error: cannot find input file '" << temp << "'. Exiting." << std::endl;
                            exit(0);
                        }
                    } else {
                        optind = curind - 1;
                        break;
                    }
                }
                logtype = "tree";
                argspresent = true;
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                argspresent = true;
                break;
            case 'b':
                burnin = string_to_int(optarg, "-b");
                if (burnin < 0) {
                    std::cerr << "Error: burnin must be a positive integer. Exiting." << std::endl;
                    exit(0);
                }
                argspresent = true;
                break;
            case 'n':
                nthin = string_to_int(optarg, "-n");
                if (nthin < 1) {
                    std::cerr << "Error: nthin must be a >= 1. Exiting." << std::endl;
                    exit(0);
                }
                argspresent = true;
                break;
            case 'r':
                nrandom = string_to_int(optarg, "-r");
                argspresent = true;
                break;
            case 'i':
                count = true;
                argspresent = true;
                break;
            case 's':
                summarize = true;
                argspresent = true;
                break;
            case 'c':
                get_columns = true;
                argspresent = true;
                break;
            case 'd':
                delete_columns = true;
                incolids = strdup(optarg);
                parse_comma_list(incolids, col_indices);
                sort(col_indices.begin(), col_indices.end());
                argspresent = true;
                break;
            case 'k':
                keep_columns = true;
                incolids = strdup(optarg);
                parse_comma_list(incolids, col_indices);
                sort(col_indices.begin(), col_indices.end());
                argspresent = true;
                break;
            case 'x':
                seed = string_to_long_int(optarg, "-x");
                argspresent = true;
                break;
            case 'v':
                verbose = true;
                argspresent = true;
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
                print_error(*argv);
                exit(0);
        }
    }
    
    if (!argspresent) {
        print_help();
        exit(1);
    }
    
    std::ostream * poos = nullptr;
    std::ofstream * ofstr = nullptr;
    
    // not used at the moment: assumed that all input comes from files
    //istream * pios = nullptr;
    //ifstream * fstr = nullptr;
    
    if (!tfileset && !pfileset) {
        std::cerr << "Error: must specify a tree file or parameter file. Exiting." << std::endl;
        exit(0);
    }
    
    if (tfileset && pfileset) {
        std::cerr << "Error: set tree file *or* parameter file, not both. Exiting." << std::endl;
        exit(0);
    }
    
    // abort if invalid args
    if (tfileset) {
        if (get_columns || delete_columns || keep_columns) {
            std::cerr << "Error: column arguments are not applicable for tree files. Exiting." << std::endl;
            exit(0);
        } else if (summarize) {
            std::cerr << "Error: summary is not applicable for tree files. Exiting." << std::endl;
            exit(0);
        }
    }
    
    // exit if not 1-indexed (just check first column
    if (delete_columns || keep_columns) {
        if (col_indices[0] < 1) {
            std::cerr << "Warning: column numbers are 1-indexed. Exiting." << std::endl;
            exit(0);
        }
    }
    
    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    //LogManipulator lm (logtype, input_files, pios, poos);
    LogManipulator lm (logtype, input_files, poos, verbose);
    
    if (count) {
        lm.count();
        lm.get_sample_counts();
    } else if (get_columns) {
        lm.get_column_names();
    } else if (delete_columns)  {
        lm.delete_columns(col_indices);
    } else if (keep_columns) {
        lm.retain_columns(col_indices);
    } else if (summarize) {
        lm.summarize(burnin, nthin);
    } else {
        lm.sample(burnin, nthin, nrandom, seed);
    }
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
