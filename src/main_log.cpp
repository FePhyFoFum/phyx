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
#include "constants.h"

extern std::string PHYX_CITATION;


void print_help() {
    std::cout << "MCMC log file manipulator." << std::endl;
    std::cout << "Can combine and resample parameters or trees across files." << std::endl;
    std::cout << "Log files need not contain the same number of samples." << std::endl;
    std::cout << "Input files may be indicated using wildcards e.g. '*.trees'" << std::endl;
    std::cout << "Parameter log files are expected to be whitespace delimited." << std::endl;
    std::cout << "*NOTE* All values are in terms of number of SAMPLES (NOT generations)." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxlog [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -p, --parmf=FILE    input parameter log file(s)" << std::endl;
    std::cout << " -t, --treef=FILE    input tree log file(s)" << std::endl;
    std::cout << " -o, --outf=FILE     output file, stout otherwise" << std::endl;
    std::cout << " -b, --burnin=INT    number of samples to exclude at the beginning of a file" << std::endl;
    std::cout << " -n, --thin=INT      interval of resampling" << std::endl;
    std::cout << " -r, --rand=INT      number of random samples (without replacement) not yet implemented!" << std::endl;
    std::cout << " -i, --info          calculate log file attributes and exit" << std::endl;
    std::cout << " -c, --columns       print out column names (parameter logs only)" << std::endl;
    std::cout << " -d, --delete=CSL    delete columns by 1-index sep by commas (NO SPACES!) (parameter logs only)" << std::endl;
    std::cout << " -k, --keep=CSL      keep only columns by 1-index sep by commas (NO SPACES!) (parameter logs only)" << std::endl;
    std::cout << " -x, --seed=INT      random number seed, clock otherwise" << std::endl;
    std::cout << " -v, --verbose       make the output more verbose" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxlog 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"parmf", required_argument, NULL, 'p'},
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"burnin", required_argument, NULL, 'b'},
    {"thin", required_argument, NULL, 'n'},
    {"rand", required_argument, NULL, 'r'},
    {"info", no_argument, NULL, 'i'},
    {"columns", no_argument, NULL, 'c'},
    {"delete", required_argument, NULL, 'd'},
    {"keep", required_argument, NULL, 'k'},
    {"seed", required_argument, NULL, 'x'},
    {"verbose", no_argument, NULL, 'v'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
    {NULL, 0, NULL, 0}
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
    int seed = -1;
    bool verbose = false;
    bool count = false;
    bool get_columns = false;
    bool delete_columns = false;
    bool keep_columns = false;
    std::string incolids;
    std::string logtype;
    
    char * outf = NULL;
    
    while (1) {
        int oi = -1;
        int curind = optind;
        int c = getopt_long(argc, argv, "p:t:o:b:n:r:icd:k:x:vhVC", long_options, &oi);
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
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'b':
                burnin = string_to_int(optarg, "-b");
                break;
            case 'n':
                nthin = string_to_int(optarg, "-n");
                break;
            case 'r':
                nrandom = string_to_int(optarg, "-r");
                break;
            case 'i':
                count = true;
                break;
            case 'c':
                get_columns = true;
                break;
            case 'd':
                delete_columns = true;
                incolids = strdup(optarg);
                parse_comma_list(incolids, col_indices);
                sort(col_indices.begin(), col_indices.end());
                break;
            case 'k':
                keep_columns = true;
                incolids = strdup(optarg);
                parse_comma_list(incolids, col_indices);
                sort(col_indices.begin(), col_indices.end());
                break;
            case 'x':
                seed = string_to_int(optarg, "-x");
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_help();
                exit(0);
            case 'V':
                std::cout << versionline << std::endl;
                exit(0);
            case 'C':
                std::cout << PHYX_CITATION << std::endl;
                exit(0);
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
    }
    
    std::ostream * poos = NULL;
    std::ofstream * ofstr = NULL;
    
    // not used at the moment: assumed that all input comes from files
    //istream * pios = NULL;
    //ifstream * fstr = NULL;
    
    if (!tfileset && !pfileset) {
        std::cerr << "Error: must specify a tree file or parameter file. Exiting." << std::endl;
        exit(0);
    }
    
    if (tfileset == true && pfileset == true) {
        std::cerr << "Error: set tree file *or* parameter file, not both. Exiting." << std::endl;
        exit(0);
    }
    
    // abort if invalid args
    if (tfileset) {
        if (get_columns || delete_columns || keep_columns) {
            std::cerr << "Error: column arguments are not applicable for tree files. Exiting." << std::endl;
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
    
    if (outfileset == true) {
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
    } else {
        lm.sample(burnin, nthin, nrandom, seed);
    }
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
