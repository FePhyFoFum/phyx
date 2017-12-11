/*
 * main_nni.cpp
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "tree_utils.h"
#include "log.h"

void print_help() {
    cout << "Nearest Neighbor Interchange Program" << endl;
    cout << "This will take newick or nexus files" << endl;
    cout << endl;
    cout << "Usage: pxnni [OPTION]... " << endl;
    cout << endl;
    cout << " -t, --treef=FILE    input tree file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output tree file, stout otherwise" << endl;
    cout << " -x, --seed=INT      random number seed, clock otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxnni 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim), Joseph F. Walker, and Joseph W. Brown");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"seed", required_argument, NULL, 'x'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    char * outf = NULL;
    char * seqf = NULL;
    int seed = -1;
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:o:x:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
            case 'x':
                seed = atoi(strdup(optarg));
                break;
            case 'V':
                cout << versionline << endl;
                exit(0);
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }

    istream * pios = NULL;
    ostream * poos = NULL;
    ifstream * fstr = NULL;
    ofstream * ofstr = NULL;
    
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    if (fileset == true) {
        fstr = new ifstream(seqf);
        pios = fstr;
    } else {
        pios = &cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }
    
    // TODO: upgrade from srand
    if (seed != -1) {
        srand(seed);
    } else {
        srand(get_clock_seed());
    }
    
    //vector<string> lines;

    string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        cerr << "this really only works with nexus or newick" << endl;
        exit(0);
    }
    
    int treeCounter = 0;
    bool going = true;
    if (ft == 1) { // newick. easy
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick (*pios, retstring, &going);
            if (tree != NULL) {
                //cout << "Working on tree #" << treeCounter << endl;
                map<Node*,vector<Node*> > tree_map;
                create_tree_map_from_rootnode(tree,tree_map);
                nni_from_tree_map(tree,tree_map);
                (*poos) << getNewickString(tree) << endl;
                delete tree;
                treeCounter++;
            }
        }
    } else if (ft == 0) { // Nexus. need to worry about possible translation tables
        map <string, string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != NULL) {
                //cout << "Working on tree #" << treeCounter << endl;
                map<Node*,vector<Node*> > tree_map;
                create_tree_map_from_rootnode(tree,tree_map);
                nni_from_tree_map(tree,tree_map);
                (*poos) << getNewickString(tree) << endl;
                delete tree;
                treeCounter++;
            }
        }
    }
    
    return EXIT_SUCCESS;
}
