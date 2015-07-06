
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <map>
using namespace std;

#include "tree.h"
#include "tree_reader.h"
#include "utils.h"

void print_help(){
    cout << "This will convert a tree file to newick." << endl;
    cout << "Can read from stdin or file." << endl;
    cout << endl;
    cout << "Usage: pxttonew [OPTION]... [FILE]..." << endl;
    cout << endl; 
    cout << " -t, --treef=FILE    input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}
/*
 * add you name if you contribute (probably add another line)
 */
string versionline("pxttonew 0.1\nCopyright (C) 2014 FePhyFoFum\nLicense GPLv2\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]){
    bool going = true;
    bool fileset = false;
    bool outfileset = false;
    char * treef;
    char * outf;
    while(going){
        int oi = -1;
        int c = getopt_long(argc,argv,"t:o:hV",long_options,&oi);
        if (c == -1){
            break;
        }
        switch(c){
            case 't':
                fileset = true;
                treef = strdup(optarg);
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
    istream * pios;
    ostream * poos;
    ifstream * fstr;
    ofstream * ofstr;
    if(fileset == true){
        fstr = new ifstream(treef);
        pios = fstr;
    }else{
        pios = &cin;
    }
    if(outfileset == true){
        ofstr = new ofstream(outf);
        poos = ofstr;
    }else{
        poos = &cout;
    }
    
    //read trees 
    string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if(ft != 0){
        cerr << "this really only converts nexus." << endl;
        exit(0);
    }
    map<string,string> translation_table;
    vector<string> retstrings;
    bool ttexists;
    ttexists = get_nexus_translation_table(*pios, &translation_table,&retstrings);
    if(retstrings.size() > 0) {
        retstring = retstrings[retstrings.size()-1];
    }
    going = true;
    Tree * tree;
    while(going){
        tree = read_next_tree_from_stream_nexus(*pios,retstring,ttexists,&translation_table, &going);
        if (going == true){
            (*poos) << tree->getRoot()->getNewick(true) << ";"<< endl;
            delete tree;
        }
    }
    if(fileset){
        fstr->close();
        delete pios;
    }if(outfileset){
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
