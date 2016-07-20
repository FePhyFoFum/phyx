
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
    cout << "This will reroot (or unroot) a tree file and produce a newick." << endl;
    cout << "Can read from stdin or file." << endl;
    cout << endl;
    cout << "Usage: pxrr [OPTION]... [FILE]..." << endl;
    cout << endl;
    cout << " -t, --treef=FILE     input tree file, stdin otherwise" << endl;
    cout << " -g, --outgroups=CSL  outgroup sep by commas (NO SPACES!)" << endl;
    cout << " -u, --unroot         unroot the tree" << endl;
    cout << " -o, --outf=FILE      output tree file, stout otherwise" << endl;
    cout << " -s, --silent         do not error if outgroup(s) not found" << endl;
    cout << " -h, --help           display this help and exit" << endl;
    cout << " -V, --version        display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxrr 0.1\nCopyright (C) 2014 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim), Joseph W. Brown");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outgroups",required_argument,NULL, 'g'},
    {"unroot", no_argument, NULL, 'u'},
    {"outf", required_argument, NULL, 'o'},
    {"silent", no_argument, NULL, 's'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outgroupsset = false;
    bool outfileset = false;
    bool silent = false;
    bool unroot = false;
    vector<string> outgroups;

    char * treef;
    char * outf;
    char * outgroupsc;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:g:uo:shV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'g':
                outgroupsset = true;
                outgroupsc = strdup(optarg);
                break;
            case 'u':
                unroot = true;
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
    if (outgroupsset == true) {
        vector<string> tokens2;
        tokenize(outgroupsc, tokens2, ",");
        for (unsigned int j=0; j < tokens2.size(); j++) {
            trim_spaces(tokens2[j]);
            outgroups.push_back(tokens2[j]);
        }
    }
    if (!outgroupsset && !unroot) {
        cerr << "you need to set the outgroup (-g)" << endl;
        exit(0);
    }

    istream * pios;
    ostream * poos;
    ifstream * fstr;
    ofstream * ofstr;
    if (fileset == true) {
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
    
    //read trees 
    string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        cerr << "this really only works with nexus or newick" << endl;
        exit(0);
    }
    bool going = true;
    bool exists;
    if (!unroot) {
        if (ft == 0) {
            map<string,string> translation_table;
            bool ttexists;
            ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);;
            Tree * tree;
            while (going) {
                tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                    &translation_table, &going);
                if (tree != NULL) {
                    exists = reroot(tree, outgroups, silent);
                    if (!exists) {
                        cerr << "the outgroup taxa don't exist in this tree " << endl;
                    } else {
                        (*poos) << tree->getRoot()->getNewick(true) << ";"<< endl;
                    }
                    delete tree;
                }
            }
        } else if (ft == 1) {
            Tree * tree;
            while (going) {
                tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
                if (going) {
                    exists = reroot(tree, outgroups, silent);
                    if (!exists) {
                        cerr << "the outgroup taxa don't exist in this tree " << endl;
                    } else {
                        (*poos) << tree->getRoot()->getNewick(true) << ";" << endl;
                    }
                    delete tree;
                }
            }
        }
    } else {
        // unroot trees
        if (ft == 0) {
            map<string,string> translation_table;
            bool ttexists;
            ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);;
            Tree * tree;
            while (going) {
                tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                    &translation_table, &going);
                if (tree != NULL) {
                    tree->unRoot();
                    (*poos) << tree->getRoot()->getNewick(true) << ";" << endl;
                    delete tree;
                }
            }
        } else if (ft == 1) {
            Tree * tree;
            while (going) {
                tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
                if (going) {
                    tree->unRoot();
                    (*poos) << tree->getRoot()->getNewick(true) << ";" << endl;
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
