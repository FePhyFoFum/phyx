#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "utils.h"
#include "log_manip.h"
#include "log.h"

void print_help() {
    cout << "MCMC log file manipulator." << endl;
    cout << "This will take parameter or tree log files." << endl;
    cout << "Can read from stdin or file." << endl;
    cout << endl;
    cout << "Usage: pxlog [OPTION]... " << endl;
    cout << endl;
    cout << " -p, --parmf=FILE    input parameter file, stdin otherwise" << endl;
    cout << " -t, --treef=FILE    input tree file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise" << endl;
    cout << " -b, --burnin=INT    number of samples to exclude at the beginning of a file" << endl;
    cout << " -n, --thin=INT      interval of resampling" << endl;
    cout << " -r, --rand=INT      number of random samples (without replacement)" << endl;
    cout << " -c, --count         just count number of samples and exit" << endl;
    cout << " -x, --seed=INT      random number seed, clock otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxlog 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"parmf", required_argument, NULL, 'p'},
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"burnin", required_argument, NULL, 'b'},
    {"thin", required_argument, NULL, 'n'},
    {"rand", required_argument, NULL, 'r'},
    {"count", required_argument, NULL, 'c'},
    {"seed", required_argument, NULL, 'x'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool pfileset = false;
    bool tfileset = false;
    
    int burnin = 0;
    int nthin = 1;
    int nrandom = -1;
    int seed = -1;
    
    bool count = false;
    
    string logtype;
    
    char * outf;
    char * treef;
    char * parmf;
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "p:t:o:b:n:r:x:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'p':
                pfileset = true;
                parmf = strdup(optarg);
                check_file_exists(parmf);
                logtype = "parm";
                break;
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
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
            case 'c':
                count = true;
                break;
            case 'x':
                seed = atoi(strdup(optarg));
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
    
    istream* pios;
    ostream* poos;
    ifstream* fstr;
    ofstream* ofstr;
    
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    if (tfileset == true && pfileset == true) {
        cout << "Set tree file *or* parameter file, not both. Exiting." << endl;
        exit (0);
    } else if (pfileset == true) {
        fstr = new ifstream(parmf);
        pios = fstr;
    } else if (tfileset == true) {
        fstr = new ifstream(treef);
        pios = fstr;
    } else {
        pios = &cin;
    }
    
    LogManipulator lm (logtype, burnin, nthin, nrandom, seed, count);
    
    /*
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        if (first) { // need to read in first sequence to get numchar
            numchar = (int)seq.get_sequence().size(); // check against partition information
            if (partitioned) {
                if (numchar != ss.get_num_partitioned_sites()) {
                    cout << "Error: numSites in sequence (" << numchar <<
                        ") does not match that in partition file (" << ss.get_num_partitioned_sites() <<
                        ")." << endl;
                }
            }
            ss.sample_sites(numchar);
            first = false;
        }
        (*poos) << ">" << seq.get_id() << endl;
        (*poos) << ss.get_resampled_seq(seq.get_sequence()) << endl;
    }
    // have to deal with last sequence outside while loop. fix this.
    // not sure if this is necessary and seems to produce duplicates at the end
    if (ft == 2) {
        (*poos) << ">" << seq.get_id() << endl;
        (*poos) << ss.get_resampled_seq(seq.get_sequence()) << endl;
    }
    */
    
    if (tfileset || pfileset) {
        fstr->close();
        delete pios;
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
