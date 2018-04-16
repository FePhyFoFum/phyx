
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <algorithm>
#include <set>
#include <map>
using namespace std;

#include "tree.h"
#include "tree_reader.h"
#include "utils.h"
#include "tree_utils.h"
#include "log.h"

void print_help() {
    cout << "This will remove tips from a tree file and produce newick" << endl;
    cout << "Can read from stdin or file" << endl;
    cout << endl;
    cout << "Usage: pxrmt [OPTION]... [FILE]..."<<endl;
    cout << endl;
    cout << " -t, --treef=FILE    input tree file, stdin otherwise" << endl;
    cout << " -n, --names=CSL     names sep by commas (NO SPACES!)" << endl;
    cout << " -f, --namesf=FILE   names in a file (each on a line)" << endl;
    cout << " -c, --comp          take the complement (i.e. move any taxa not in list)" << endl;
    cout << " -o, --outf=FILE     output tree file, stout otherwise" << endl;
    cout << " -s, --silent        suppress warnings of missing tips" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "NOTE: if you get a segfault, you may try unrooting (pxrr -u) before pruning" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" <<endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>"<<endl;
}
/*
 * add you name if you contribute (probably add another line)
 */
string versionline("pxrmt 0.1\nCopyright (C) 2014 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim), Joseph W. Brown");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"names",required_argument,NULL,'n'},
    {"outf", required_argument, NULL, 'o'},
    {"comp", no_argument, NULL, 'c'},
    {"silent", required_argument, NULL, 's'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
      
    bool fileset = false;
    bool namesset = false;
    bool namefileset = false;
    bool outfileset = false;
    bool silent = false;
    vector<string> names;
    bool complement = false;

    char * treef = NULL;
    char * outf = NULL;
    char * namesc = NULL;
    char * namesfc = NULL;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:n:f:co:shV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'n':
                namesset = true;
                namesc = strdup(optarg);
                break;
            case 'f':
                namefileset = true;
                namesfc = strdup(optarg);
                check_file_exists(namesfc);
                break;
            case 'c':
                complement = true;
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 's':
                silent = true;
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
    
    if (namesset == true) {
        vector<string> tokens2;
        string del2(",");
        tokens2.clear();
        tokenize(namesc, tokens2, del2);
        for (unsigned int j=0; j < tokens2.size(); j++) {
            trim_spaces(tokens2[j]);
            names.push_back(tokens2[j]);
        }
    } else if (namefileset == true) {
        ifstream nfstr(namesfc);
        string tline;
        while (getline(nfstr,tline)) {
            trim_spaces(tline);
            names.push_back(tline);
        }
        nfstr.close();
    } else {
        cerr << "you need to set the names of the tips you want to remove (-n)" << endl;
        exit(0);
    }

    istream * pios = NULL;
    ostream * poos = NULL;
    ifstream * fstr = NULL;
    ofstream * ofstr = NULL;
    
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
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    //read trees 
    string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        cerr << "this really only works with nexus or newick" << endl;
        exit(0);
    }
    
    bool going = true;
    if (!complement) {
        if (ft == 0) {
            map<string,string> translation_table;
            bool ttexists;
            ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
            Tree * tree;
            while (going) {
                tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                    &translation_table, &going);
                if (going == true) {
                    remove_tips(tree, names, silent);
                    (*poos) << getNewickString(tree) << endl;
                    delete tree;
                }
            }
        } else if (ft == 1) {
            Tree * tree;
            while (going) {
                tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
                if (going == true) {
                    remove_tips(tree, names, silent);
                    (*poos) << getNewickString(tree) << endl;
                    delete tree;
                }
            }
        }
    } else {
        // *** check list of names to keep is at least 2
        if (ft == 0) {
            map<string,string> translation_table;
            bool ttexists;
            ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
            Tree * tree;
            int numLeaves;
            vector <string> toPrune;
            while (going) {
                tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                    &translation_table, &going);
                if (going == true) {
                    toPrune = get_complement_tip_set(tree, names);
                    numLeaves = tree->getExternalNodeCount();
                    if (numLeaves - (int)toPrune.size() > 1) {
                        remove_tips(tree, toPrune, silent);
                        (*poos) << getNewickString(tree) << endl;
                    }
                    delete tree;
                }
            }
        } else if (ft == 1) {
            Tree * tree;
            vector <string> toPrune;
            int numLeaves;
            while (going) {
                tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
                if (going == true) {
                    toPrune = get_complement_tip_set(tree, names);
                    numLeaves = tree->getExternalNodeCount();
                    //cout << "numLeaves = " << numLeaves << endl;
                    //print_vector(toPrune);
                    if (numLeaves - (int)toPrune.size() > 1) {
                        remove_tips(tree, toPrune, silent);
                        (*poos) << getNewickString(tree) << endl;
                    }
                    delete tree;
                }
            }
        }
    }
    
    if (fileset) {
        fstr->close();
        delete pios;
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
