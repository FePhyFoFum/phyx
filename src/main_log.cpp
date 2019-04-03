#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "utils.h"
#include "log_manip.h"
#include "log.h"

void print_help() {
    cout << "MCMC log file manipulator." << endl;
    cout << "Can combine and resample parameters or trees across files." << endl;
    cout << "Log files need not contain the same number of samples." << endl;
    cout << "Input files may be indicated using wildcards e.g. '*.trees'" << endl;
    cout << "Parameter log files are expected to be whitespace delimited." << endl;
    cout << "*NOTE* All values are in terms of number of SAMPLES (NOT generations)." << endl;
    cout << endl;
    cout << "Usage: pxlog [OPTION]... " << endl;
    cout << endl;
    cout << " -p, --parmf=FILE    input parameter log file(s)" << endl;
    cout << " -t, --treef=FILE    input tree log file(s)" << endl;
    cout << " -o, --outf=FILE     output file, stout otherwise" << endl;
    cout << " -b, --burnin=INT    number of samples to exclude at the beginning of a file" << endl;
    cout << " -n, --thin=INT      interval of resampling" << endl;
    cout << " -r, --rand=INT      number of random samples (without replacement) not yet implemented!" << endl;
    cout << " -i, --info          calculate log file attributes and exit" << endl;
    cout << " -c, --columns       print out column names (parameter logs only)" << endl;
    cout << " -d, --delete=CSL    delete columns by 1-index sep by commas (NO SPACES!) (parameter logs only)" << endl;
    cout << " -k, --keep=CSL      keep only columns by 1-index sep by commas (NO SPACES!) (parameter logs only)" << endl;
    cout << " -x, --seed=INT      random number seed, clock otherwise" << endl;
    cout << " -v, --verbose       make the output more verbose" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxlog 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

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
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool pfileset = false;
    bool tfileset = false;
    vector <string> input_files;
    vector <int> col_indices;
    
    int burnin = 0;
    int nthin = 1;
    int nrandom = -1;
    int seed = -1;
    bool verbose = false;
    bool count = false;
    bool get_columns = false;
    bool delete_columns = false;
    bool keep_columns = false;
    string incolids;
    
    string logtype;
    
    char * outf = NULL;
    
    while (1) {
        int oi = -1;
        int curind = optind;
        int c = getopt_long(argc, argv, "p:t:o:b:n:r:icd:k:x:vhV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'p':
                pfileset = true;
                curind = optind - 1;
                while (curind < argc) {
                    string temp = strdup(argv[curind]);
                    curind++;
                    if (temp[0] != '-') {
                        ifstream infile(temp.c_str());
                        if (infile.good()) { // check that file exists
                            input_files.push_back(temp);
                            infile.close();
                        } else {
                            cout << "Cannot find input file '" << temp << "'. Exiting." << endl;
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
                    string temp = strdup(argv[curind]);
                    curind++;
                    if (temp[0] != '-') {
                        ifstream infile(temp.c_str());
                        if (infile.good()) { // check that file exists
                            input_files.push_back(temp);
                            infile.close();
                        } else {
                            cout << "Cannot find input file '" << temp << "'. Exiting." << endl;
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
                burnin = atoi(strdup(optarg));
                break;
            case 'n':
                nthin = atoi(strdup(optarg));
                break;
            case 'r':
                nrandom = atoi(strdup(optarg));
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
                seed = atoi(strdup(optarg));
                break;
            case 'v':
                verbose = true;
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
    
    ostream * poos = NULL;
    ofstream * ofstr = NULL;
    
    // not used at the moment: assumed that all input comes from files
    //istream * pios = NULL;
    //ifstream * fstr = NULL;
    
    if (!tfileset && !pfileset) {
        cout << "Must specify a tree file or parameter file. Exiting." << endl;
        exit (0);
    }
    
    if (tfileset == true && pfileset == true) {
        cout << "Set tree file *or* parameter file, not both. Exiting." << endl;
        exit (0);
    }
    
    // abort if invalid args
    if (tfileset) {
        if (get_columns || delete_columns || keep_columns) {
            cout << "Column arguments are not applicable for tree files. Exiting." << endl;
            exit (0);
        }
    }
    
    // exit if not 1-indexed (just check first column
    if (delete_columns || keep_columns) {
        if (col_indices[0] < 1) {
            cout << "Warning: column numbers are 1-indexed. Exiting." << endl;
            exit (0);
        }
    }
    
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
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
