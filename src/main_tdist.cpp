// not currently used

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

// altering bd code. doesn't currently do anything!

/*
Give two options:
1. pass in 2 trees
2. pass in 1 distribution of trees
 */

void print_help () {
    std::cout << "Calculate tree distances. RF to begin with, others to follow" << std::endl;
    std::cout << "Either pass in 2 trees with `t` and `a`, or a single distribution with `t`." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxtdist [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -t, --treef=FILE    reference treefile, stdin otherwise" << std::endl;
    std::cout << " -a, --alttree=FILE  alternate treefile" << std::endl;
    std::cout << " -d, --dist=STRING   distance metric, default='RF'" << std::endl;
    std::cout << " -o, --outf=FILE     output file, stout otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

string versionline("pxtdist 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"alttree", required_argument, NULL, 'a'},
    {"dist", required_argument, NULL, 'd'},
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
    bool alttfileset = false;
    
    char * treef = NULL;
    char * alttreef = NULL;
    char * outf = NULL;
    
    std::string dist = "bd";
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:a:d:o:x:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'a':
                alttfileset = true;
                alttreef = strdup(optarg);
                check_file_exists(alttreef);
                break;
            case 'd':
                dist = strdup(optarg);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
            case 'V':
                std::cout << versionline << std::endl;
                exit(0);
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
    }
    
    if (tfileset && outfileset) {
        check_inout_streams_identical(treef, outf);
    }
    if (alttfileset && outfileset) {
        check_inout_streams_identical(alttreef, outf);
    }
    
    std::istream * pios = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ofstream * ofstr = NULL;

    if (outfileset == true) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    if (tfileset == true) {
        fstr = new std::ifstream(treef);
        pios = fstr;
    } else {
        pios = &cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }
    
    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        cerr << "this really only works with nexus or newick" << std::endl;
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
                std::cout << "Tree is not ultrametric. Exiting." << std::endl;
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
