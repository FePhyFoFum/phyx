
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "utils.h"
#include "tree_reader.h"
#include "rename.h"
#include "log.h"

void print_help() {
    cout << "Taxon relabelling for trees." << endl;
    cout << "This will take nexus and newick inputs." << endl;
    cout << "Two ordered lists of taxa, -c (current) and -n (new) must be provided." << endl;
    cout << endl;
    cout << "Usage: pxrnt [OPTION]... " << endl;
    cout << endl;
    cout << " -t, --treef=FILE    input tree file, stdin otherwise" << endl;
    cout << " -c, --cnames=FILE   file containing current taxon labels (one per line)" << endl;
    cout << " -n, --nnames=FILE   file containing new taxon labels (one per line)" << endl;
    cout << " -o, --outf=FILE     output file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxrnt 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 'f'},
    {"cnames", required_argument, NULL, 'c'},
    {"nnames", required_argument, NULL, 'n'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool tfileset = false;
    bool cfileset = false;
    bool nfileset = false;
    char * outf = NULL;
    char * treef = NULL;
    string cnamef = "";
    string nnamef = "";
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:c:n:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'c':
                cfileset = true;
                cnamef = strdup(optarg);
                check_file_exists(cnamef.c_str());
                break;
            case 'n':
                nfileset = true;
                nnamef = strdup(optarg);
                check_file_exists(nnamef.c_str());
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
    
    if (!nfileset | !cfileset) {
        cout << "Must supply both name files (-c for current, -n for new)." << endl;
        exit(0);
    }
    
    if (tfileset == true) {
        fstr = new ifstream(treef);
        pios = fstr;
    } else {
        pios = &cin;
    }
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    Rename rn (cnamef, nnamef);
    
    string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        cerr << "this really only works with nexus or newick" << endl;
        exit(0);
    }
    bool going = true;
    if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (going) {
                rn.rename_tree(tree);
                //exists = reroot(tree, outgroups, silent);
                //if (!exists) {
                //    cerr << "the outgroup taxa don't exist in this tree " << endl;
                //} else {
                    (*poos) << tree->getRoot()->getNewick(true) << ";" << endl;
                //}
                delete tree;
            }
        }
    }
    
    /*
    if (fileset == true) {
        fstr = new ifstream(seqf);
        pios = fstr;
    } else {
        pios = &cin;
    }
    
    SequenceRecoder sr;
    
    Sequence seq;
    string retstring;
    //bool first = true;
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        (*poos) << ">" << seq.get_id() << endl;
        (*poos) << sr.get_recoded_seq(seq.get_sequence()) << endl;
    }
// have to deal with last sequence outside while loop. fix this.
    if (ft == 2) {
        (*poos) << ">" << seq.get_id() << endl;
        (*poos) << sr.get_recoded_seq(seq.get_sequence()) << endl;
    }
    if (fileset) {
        fstr->close();
        delete pios;
    }
    */
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
