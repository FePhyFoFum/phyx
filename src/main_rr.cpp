
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
    cout << " -r, --ranked         turn on ordering of outgroups. will root on first one present" << endl;
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
    {"outgroups", required_argument, NULL, 'g'},
    {"ranked", no_argument, NULL, 'r'},
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
    bool ranked = false;
    vector <string> outgroups;

    char * treef = NULL;
    char * outf = NULL;
    char * outgroupsc = NULL;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:g:ruo:shV", long_options, &oi);
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
            case 'r':
                ranked = true;
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
    
    if (fileset && outfileset) {
        check_inout_streams_identical(treef, outf);
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
    bool exists;
    if (!unroot) {
        if (ft == 0) { // Nexus
            map<string,string> translation_table;
            bool ttexists;
            ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);;
            Tree * tree;
            while (going) {
                tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                    &translation_table, &going);
                if (tree != NULL) {
                    if (ranked) {
                        // find first outgroup present in tree
                        bool ogexists = false;
                        for (unsigned int i=0; i < outgroups.size(); i++) {
                            string name = outgroups[i];
                            if (check_name_against_tree(tree, name)) {
                                vector <string> og;
                                og.push_back(name);
                                exists = reroot(tree, og, silent);
                                ogexists = true;
                                break;
                            }
                        }
                        // if no valid outgroups, let silent option figure out
                        if (!ogexists) {
                            exists = reroot(tree, outgroups, silent);
                        }
                    } else {
                        exists = reroot(tree, outgroups, silent);
                    }
                    if (!exists) {
                        cerr << "the outgroup taxa don't exist in this tree " << endl;
                    } else {
                        (*poos) << getNewickString(tree) << endl;
                    }
                    delete tree;
                }
            }
        } else if (ft == 1) { // newick
            Tree * tree;
            while (going) {
                tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
                if (going) {
                    if (ranked) {
                        // find first outgroup present in tree
                        bool ogexists = false;
                        for (unsigned int i=0; i < outgroups.size(); i++) {
                            string name = outgroups[i];
                            if (check_name_against_tree(tree, name)) {
                                vector <string> og;
                                og.push_back(name);
                                exists = reroot(tree, og, silent);
                                ogexists = true;
                                break;
                            }
                        }
                        // if no valid outgroups, let silent option figure out
                        if (!ogexists) {
                            exists = reroot(tree, outgroups, silent);
                        }
                    } else {
                        exists = reroot(tree, outgroups, silent);
                    }
                    if (!exists) {
                        cerr << "the outgroup taxa don't exist in this tree " << endl;
                    } else {
                        (*poos) << getNewickString(tree) << endl;
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
                    (*poos) << getNewickString(tree) << endl;
                    delete tree;
                }
            }
        } else if (ft == 1) {
            Tree * tree;
            while (going) {
                tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
                if (going) {
                    tree->unRoot();
                    (*poos) << getNewickString(tree) << endl;
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
