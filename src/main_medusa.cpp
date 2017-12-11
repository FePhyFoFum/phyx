#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>
#include <algorithm>
#include <cmath>

using namespace std;

#include "tree_reader.h"
#include "tree.h"
#include "tree_utils.h"
#include "utils.h"
#include "bd_fit.h"
#include "log.h"

void print_help () {
    cout << "Piecewise birth-death inference using AIC" << endl;
    cout << endl;
    cout << "Usage: pxbdfit [OPTION]... " << endl;
    cout << endl;
    cout << " -t, --treef=FILE    input treefile, stdin otherwise" << endl;
    cout << " -m, --model=STRING  diversification model; either 'yule' or 'bd' (default)" << endl;
    cout << " -o, --outf=FILE     output file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxbdfit 0.1\nCopyright (C) 2017 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"model", required_argument, NULL, 'm'},
    {"outf", required_argument, NULL, 'o'},
    {"showd", no_argument, NULL, 's'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool tfileset = false;
    
    char * treef = NULL;
    char * outf = NULL;
    
    string model = "bd";
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:m:o:x:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'm':
                // need to check valid models here
                model = strdup(optarg);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
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
    
    istream* pios = NULL;
    ostream* poos = NULL;
    ifstream* fstr = NULL;
    ofstream* ofstr = NULL;

    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    if (tfileset == true) {
        fstr = new ifstream(treef);
        pios = fstr;
    } else {
        pios = &cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }
    
    string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        cerr << "this really only works with nexus or newick" << endl;
        exit(0);
    }
    
    bool going = true;
    Tree * tree;
    while (going) {
        tree = read_next_tree_from_stream_newick (*pios, retstring, &going);
        if (going) {
            // in addition to checking ultramtericity, the following sets node heights
            //if (is_ultrametric_postorder(tree)) {
            if (is_ultrametric_paths(tree)) {
                BDFit bd(tree, model);
                bd.get_pars(poos);
                delete tree;
            } else {
                cout << "Tree is not ultrametric. Exiting." << endl;
                exit(0);
            }
        }
    }
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
