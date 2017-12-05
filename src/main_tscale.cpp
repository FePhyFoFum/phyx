
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "utils.h"
#include "tree_reader.h"
#include "tscale.h"
#include "tree_utils.h"
#include "log.h"

void print_help() {
    cout << "Tree rescaling by providing either scaling factor or root height." << endl;
    cout << "This will take nexus and newick inputs." << endl;
    cout << endl;
    cout << "Usage: pxtscale [OPTION]... " << endl;
    cout << endl;
    cout << " -t, --treef=FILE        input tree file, stdin otherwise" << endl;
    cout << " -s, --scale=DOUBLE      edge length scaling factor" << endl;
    cout << " -r, --rootheight=DOUBLE height of root (tree must be ultrametric)" << endl;
    cout << " -o, --outf=FILE         output file, stout otherwise" << endl;
    cout << " -h, --help              display this help and exit" << endl;
    cout << " -V, --version           display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxtscale 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"scale", required_argument, NULL, 's'},
    {"rootheight", required_argument, NULL, 'r'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool tfileset = false;
    double rootheight = 0.0;
    double scalef = 1.0;
    bool heightset = false;
    bool scaleset = false;
    char * outf = NULL;
    char * treef = NULL;
    string cnamef = "";
    string nnamef = "";
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:s:r:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 's':
                scaleset = true;
                scalef = atof(strdup(optarg));
                break;
            case 'r':
                heightset = true;
                rootheight = atof(strdup(optarg));
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
    
    if (tfileset && outfileset) {
        check_inout_streams_identical(treef, outf);
    }
    
    istream * pios = NULL;
    ostream * poos = NULL;
    ifstream * fstr = NULL;
    ofstream * ofstr = NULL;
    
    if (heightset && scaleset) {
        cout << "Supply only 'rootheight' or 'scale', not both. Exiting." << endl;
        exit(0);
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
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    TScale ts;
    
    if (heightset) {
        ts.set_rootheight(rootheight);
    } else {
        ts.set_scalef(scalef);
    }
    
    string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        cerr << "This really only works with nexus or newick" << endl;
        exit(0);
    }
    bool going = true;
    if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (going) {
                if (heightset) {
                    // have to check ultrametricity
                    bool isultra = is_ultrametric_paths(tree);
                    if (!isultra) {
                        cout << "Setting root height only works for ultrametric trees. Exiting." << endl;
                        exit(0);
                    }
                }
                ts.rescale(tree);
                (*poos) << getNewickString(tree) << endl;
                delete tree;
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
                if (heightset) {
                    // have to check ultrametricity
                    bool isultra = is_ultrametric_paths(tree);
                    if (!isultra) {
                        cout << "Setting root height only works for ultrametric trees. Exiting." << endl;
                        exit(0);
                    }
                }
                ts.rescale(tree);
                (*poos) << getNewickString(tree) << endl;
                delete tree;
            }
        }
    }
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
