
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

void print_help() {
    cout << "This will remove tips from a tree file and produce newick" << endl;
    cout << "Can read from stdin or file" << endl;
    cout << endl;
    cout << "Usage: pxrmt [OPTION]... [FILE]..."<<endl;
    cout << endl; 
    cout << " -t, --treef=FILE     input tree file, stdin otherwise"<<endl;
    cout << " -n, --names=CSL      names sep by commas (NO SPACES!)" <<endl;
    cout << " -f, --namesf=FILE    names in a file (each on a line)" << endl;
    cout << " -o, --outf=FILE      output tree file, stout otherwise"<<endl;
    cout << "     --help           display this help and exit"<<endl;
    cout << "     --version        display version and exit"<<endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" <<endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>"<<endl;
}
/*
 * add you name if you contribute (probably add another line)
 */
string versionline("pxrmt 0.1\nCopyright (C) 2014 FePhyFoFum\nLicense GPLv2\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"names",required_argument,NULL,'n'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

bool removetip(Tree * tree, vector<string> & names);
bool removetip(Tree * tree,vector<string> & names) {
    for (int i=0; i < names.size(); i++) {
        Node * m = tree->getExternalNode(names[i]);
        if (m != NULL) {
            tree->pruneExternalNode(m);
        } else {
            cerr << names[i] << " not in tree"  << endl;
        }
    }
    return true;
}

int main(int argc, char * argv[]) {
    bool fileset = false;
    bool namesset = false;
    bool namefileset = false;
    bool outfileset = false;
    vector<string> names;

    char * treef;
    char * outf;
    char * namesc;
    char * namesfc;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:n:f:o:hV", long_options, &oi);
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
    if (ft == 0) {
        map<string,string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (going == true) {
                exists = removetip(tree,names);
                (*poos) << tree->getRoot()->getNewick(true) << ";" << endl;
                delete tree;
            }
        }
    } else if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (going == true) {
                exists = removetip(tree, names);
                if (exists == false) {
                    cerr << "the names don't exist in this tree " << endl;
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
