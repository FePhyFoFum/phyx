/*
 * main_bd.cpp
 *
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "bd_sim.h"
#include "log.h"

void print_help () {
    cout << "Birth-death simulator." << endl;
    cout << endl;
    cout << "Usage: pxbdsim [OPTION]... " << endl;
    cout << endl;
    cout << " -e, --extant=INT    number of extant species, alt to time" << endl;
    cout << " -t, --time=INT      depth of the tree, alt to extant" << endl;
    cout << " -b, --birth=DOUBLE  birth rate, default=1" << endl;
    cout << " -d, --death=DOUBLE  death rate, default=0" << endl;
    cout << " -n, --nreps=INT     number of replicates, default=1" << endl;
    cout << " -o, --outf=FILE     output file, stout otherwise" << endl;
    cout << " -s, --showd         show dead taxa" << endl;
    cout << " -x, --seed=INT      random number seed, clock otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxbdsim 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim), Joseph W. Brown");

static struct option const long_options[] =
{
    {"extant", required_argument, NULL, 'e'},
    {"time", required_argument, NULL, 't'},
    {"birth", required_argument, NULL, 'b'},
    {"death", required_argument, NULL, 'd'},
    {"nreps", required_argument, NULL, 'n'},
    {"outf", required_argument, NULL, 'o'},
    {"showd", no_argument, NULL, 's'},
    {"seed", required_argument, NULL, 'x'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool timeset = false;
    bool extantset = false;
    char * outf = NULL;
    int ext = 0;
    int nreps = 1;
    double time = 0.0;
    double birth = 1.0;
    double death = 0.0;
    bool showd = false;
    
    int seed = -1;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "e:t:b:d:n:o:x:shV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'e':
                ext = atoi(strdup(optarg));
                extantset = true;
                break;
            case 't':
                time = atof(strdup(optarg));
                timeset = true;
                break;
            case 'b':
                birth = atof(strdup(optarg));
                if (birth <= 0) {
                    cout << "Birth rate must be > 0" << endl;
                    exit(0);
                }
                break;
            case 'd':
                death = atof(strdup(optarg));
                if (death < 0) {
                    cout << "Death rate must be >= 0" << endl;
                    exit(0);
                }
                break;
            case 'n':
                nreps = atoi(strdup(optarg));
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'x':
                seed = atoi(strdup(optarg));
                break;
            case 's':
                showd = true;
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
    
    if (ext == 0 && time == 0) {
        cout << "you have to set -e or -t" << endl;
        exit(0);
    }
    if (timeset && extantset) {
        cout << "Set -e or -t, not both" << endl;
        exit(0);
    }
    
    ostream * poos = NULL;
    ofstream * ofstr = NULL;
    
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    BirthDeathSimulator bd(ext, time, birth, death, seed);
    for (int i = 0; i < nreps; i++) {
        Tree * bdtr = bd.make_tree(showd);
        if (bdtr->getExtantNodeCount() > 1) {
            (*poos) << bdtr->getRoot()->getNewick(true) << ";" << endl;
        } else {
            (*poos) << "(" << bdtr->getRoot()->getNewick(true) << ");" << endl;
        }
        
        
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
