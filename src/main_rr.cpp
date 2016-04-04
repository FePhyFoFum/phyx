
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

void print_help(){
    cout << "This will reroot a tree file and produce a newick." << endl;
    cout << "Can read from stdin or file." << endl;
    cout << endl;
    cout << "Usage: pxrr [OPTION]... [FILE]..." << endl;
    cout << endl; 
    cout << " -t, --treef=FILE     input tree file, stdin otherwise" << endl;
    cout << " -g, --outgroups=CSL  outgroup sep by commas (NO SPACES!)" << endl; 
    cout << " -o, --outf=FILE      output tree file, stout otherwise" << endl;
    cout << "     --help           display this help and exit" << endl;
    cout << "     --version        display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}
/*
 * add you name if you contribute (probably add another line)
 */
string versionline("pxrr 0.1\nCopyright (C) 2014 FePhyFoFum\nLicense GPLv2\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outgroups",required_argument,NULL, 'g'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

bool reroot(Tree * tree, vector<string> & outgr);
bool reroot(Tree * tree, vector<string> & outgr) {
    if (!check_names_against_tree(tree, outgr)) {
        return false;
    }
    Node * m = tree->getMRCA(outgr);
    if (m == NULL) {
        return false;
    }
    if (m == tree->getRoot()) {
        //check to see if the outgroups are just the children of the root
        //if so, then do this
        //tree->rootWithRootTips(outgr);
        //if not, then do this
        tree->getInternalMRCA(outgr);
        return true;
    }
    bool success = tree->reRoot(m);
    return success;
}

int main(int argc, char * argv[]) {
    bool going = true;
    bool fileset = false;
    bool outgroupsset = false;
    bool outfileset = false;
    vector<string> outgroups;

    char * treef;
    char * outf;
    char * outgroupsc;
    while (going) {
        int oi = -1;
        int c = getopt_long(argc, argv, "g:t:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                break;
            case 'g':
                outgroupsset = true;
                outgroupsc = strdup(optarg);
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
                print_error(argv[0],(char)c);
                exit(0);
        }
    }
    if (outgroupsset == true) {
        vector<string> tokens2;
        string del2(",");
        tokens2.clear();
        tokenize(outgroupsc, tokens2, del2);
        for (unsigned int j=0; j < tokens2.size(); j++) {
            trim_spaces(tokens2[j]);
            outgroups.push_back(tokens2[j]);
        }
    } else {
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
    going = true;
    bool exists;
    if (ft == 0) {
        map<string,string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);;
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != NULL) {
                exists = reroot(tree, outgroups);
                (*poos) << tree->getRoot()->getNewick(true) << ";"<< endl;
                delete tree;
            }
        }
    } else if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (going == true) {
                exists = reroot(tree, outgroups);
                if (exists == false) {
                    cerr << "the outgroup taxa don't exist in this tree " << endl;
                } else {
                    (*poos) << tree->getRoot()->getNewick(true) << ";" << endl;
                }
                delete tree;
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
