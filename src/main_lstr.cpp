#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "tree_info.h"
#include "tree_reader.h"
#include "tree.h"
#include "utils.h"
#include "log.h"

void print_help() {
    cout << "Print tree summary" << endl;
    cout << "By default returns all properties. Alternatively choose 1 property." << endl;
    cout << "This will take newick or nexus files" << endl;
    cout << endl;
    cout << "Usage: pxlstr [OPTION]... " << endl;
    cout << endl;
    cout << " -t, --treef=FILE    input tree file, stdin otherwise" << endl;
    cout << " -r, --rooted        return whether the tree is rooted" << endl;
    cout << " -a, --age           return the height of root (must be rooted and ultrametric)" << endl;
    cout << " -n, --ntips         return the number of terminals" << endl;
    cout << " -u, --ultrametric   return whether tree is ultrametric" << endl;
    cout << " -b, --binary        return whether tree is binary" << endl;
    cout << " -l, --length        return the length of the tree" << endl;
    cout << " -i, --tiplabels     return all tip labels (one per line)" << endl;
    cout << " -v, --rtvar         return root-to-tip variance" << endl;
    cout << " -o, --outf=FILE     output tree stats file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxlstr 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"rooted", no_argument, NULL, 'r'},
    {"age", no_argument, NULL, 'a'},
    {"ntips", no_argument, NULL, 'n'},
    {"ultrametric", no_argument, NULL, 'u'},
    {"binary", no_argument, NULL, 'b'},
    {"length", no_argument, NULL, 'l'},
    {"tiplabels", no_argument, NULL, 'i'},
    {"rtvar",no_argument, NULL, 'v'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    bool optionsset = false; // if true, do not return all properties
    bool ultracheck = false;
    bool binarycheck = false;
    bool lengthcheck = false;
    bool agecheck = false;
    bool rootedcheck = false;
    bool ntipcheck = false;
    bool namecheck = false;
    bool rtvarcheck = false;
    char * outf = NULL;
    char * treef = NULL;
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:vranublio:x:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'r':
                rootedcheck = true;
                optionsset = true;
                break;
            case 'a':
                agecheck = true;
                optionsset = true;
                break;
            case 'n':
                ntipcheck = true;
                optionsset = true;
                break;
            case 'u':
                ultracheck = true;
                optionsset = true;
                break;
            case 'b':
                binarycheck = true;
                optionsset = true;
                break;
            case 'l':
                lengthcheck = true;
                optionsset = true;
                break;
            case 'i':
                namecheck = true;
                optionsset = true;
                break;
            case 'v':
                rtvarcheck = true;
                optionsset = true;
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
    
    if (fileset && outfileset) {
        check_inout_streams_identical(treef, outf);
    }
    
    if ((ultracheck + binarycheck + lengthcheck + agecheck + rootedcheck + rtvarcheck + ntipcheck + namecheck) > 1) {
        cout << "Specify 1 property only (or leave blank to show all properties)" << endl;
        exit(0);
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
    
    int treeCounter = 0;
    bool going = true;
    if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (tree != NULL) {
                if (!optionsset) {
                    (*poos) << "tree #: " << treeCounter << endl;
                    TreeInfo ti(tree);
                    ti.get_stats(poos);
                    delete tree;
                    treeCounter++;
                } else {
                    // only a single property
                    TreeInfo ti(tree, ultracheck, binarycheck, agecheck, rootedcheck,
                        ntipcheck, lengthcheck, namecheck, rtvarcheck, poos);
                }
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
                if (!optionsset) {
                    (*poos) << "tree #: " << treeCounter << endl;
                    TreeInfo ti(tree);
                    ti.get_stats(poos);
                    delete tree;
                    treeCounter++;
                } else {
                    // only a single property
                    TreeInfo ti(tree, ultracheck, binarycheck, agecheck, rootedcheck,
                        ntipcheck, lengthcheck, namecheck, rtvarcheck, poos);
                }
            }
        }
    }
    
    return EXIT_SUCCESS;
}
