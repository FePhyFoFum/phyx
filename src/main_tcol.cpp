
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <map>
#include <getopt.h>

using namespace std;

#include "utils.h"
#include "tree_reader.h"
#include "tree_utils.h"
#include "log.h"

void print_help() {
    cout << "Add information to a tree so that you can color the edges." << endl;
    cout << "This will take nexus and newick inputs and will output nexus so" << endl;
    cout << "that it can be read by figtree." << endl;
    cout << endl;
    cout << "Usage: pxtcol [OPTION]... " << endl;
    cout << endl;
    cout << " -t, --treef=FILE        input tree file, stdin otherwise" << endl;
    cout << " -m, --mrcaf=FILE        file with mrcas and annotations, tab separated" << endl;
    cout << " -d, --nodeidf=FILE      file with nodeids (labels) and annotations, tab separated" << endl;
    cout << " -o, --outf=FILE         output file, stout otherwise" << endl;
    cout << " -h, --help              display this help and exit" << endl;
    cout << " -V, --version           display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxtcol 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"mrcaf", required_argument, NULL, 's'},
    {"nodeidf", required_argument, NULL, 'r'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool tfileset = false;
    bool mrcafset = false;
    bool nodeidfset = false;
    char * mrcaf = NULL;
    char * nodeidf = NULL;
    char * outf = NULL;
    char * treef = NULL;
    string cnamef = "";
    string nnamef = "";
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:m:d:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'm':
                mrcafset = true;
                mrcaf = strdup(optarg);
                break;
            case 'd':
                nodeidfset = true;
                nodeidf = strdup(optarg);
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
    
    //for node ids
    map<string, string> nodeid_map;
    if (nodeidfset == true) {
        ifstream nfstr(nodeidf);
        string tline;
        while (getline(nfstr,tline)) {
            trim_spaces(tline);
            vector<string> tokens2;
            tokenize(tline, tokens2, "\t");
            for (unsigned int j=0; j < tokens2.size(); j++) {
                trim_spaces(tokens2[j]);
            }
            if (tokens2.size() != 2)
                continue;
            nodeid_map[tokens2[0]] = tokens2[1];
        }
        nfstr.close();
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
        (*poos) << "#NEXUS" << endl;
        (*poos) << "begin trees;" << endl;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (going) {
                for (int i=0; i < tree->getInternalNodeCount(); i++) {
                    Node * tnode = tree->getInternalNode(i);
                    if (nodeid_map.find(tnode->getName()) != nodeid_map.end()) {
                        tnode->setName("[&name=\""+tnode->getName()+"\",ann="+nodeid_map[tnode->getName()]+"]");
                    } 
                }
                for (int i=0; i < tree->getExternalNodeCount(); i++) {
                    Node * tnode = tree->getExternalNode(i);
                    if (nodeid_map.find(tnode->getName()) != nodeid_map.end()) {
                        tnode->setName(tnode->getName()+"[&ann="+nodeid_map[tnode->getName()]+"]");
                    } 
                }
                // put annotations here
                (*poos) << "tree tree = " << getNewickString(tree) << endl;
                delete tree;
            }
        }
        (*poos) << "end;" << endl;
    } else if (ft == 0) { // Nexus. need to worry about possible translation tables
        map <string, string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        Tree * tree;
        (*poos) << "#NEXUS" << endl;
        (*poos) << "begin trees;" << endl;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != NULL) {
                for (int i=0; i < tree->getInternalNodeCount(); i++) {
                    Node * tnode = tree->getInternalNode(i);
                    if (nodeid_map.find(tnode->getName()) != nodeid_map.end()){
                        tnode->setName("[&name=\""+tnode->getName()+"\",ann="+nodeid_map[tnode->getName()]+"]");
                    } 
                }
                for (int i=0; i < tree->getExternalNodeCount(); i++) {
                    Node * tnode = tree->getExternalNode(i);
                    if (nodeid_map.find(tnode->getName()) != nodeid_map.end()){
                        tnode->setName(tnode->getName()+"[&ann="+nodeid_map[tnode->getName()]+"]");
                    } 
                }   
                // put annotations here
                (*poos) << "tree tree = " << getNewickString(tree) << endl;
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
