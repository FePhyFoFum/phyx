/*
 * The idea behind this is to allow for the naming of internal nodes based
 * on given MRCAS and a set of names in a file the input of which should
 * look like:
 * MRCANAME = tip1 tip2 ...
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

#include "node.h"
#include "tree_reader.h"
#include "tree.h"
#include "utils.h"
#include "tree_utils.h"
#include "log.h"

void print_help() {
    std::cout << "Label internal nodes with clade names." << std::endl;
    std::cout << "Takes in newick tree and MRCA file with format:" << std::endl;
    std::cout << "MRCANAME = tip1 tip2 ..." << std::endl;
    std::cout << "If no MRCA file is present, this will label anything" << std::endl;
    std::cout << "that isn't labeled" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxmrcaname [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -t, --treef=FILE    input newick tree file, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output newick file, stout otherwise" << std::endl;
    std::cout << " -m, --mrca=FILE     file containing MRCA declarations" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxmrcaname 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim), Joseph W. Brown");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"mrca", required_argument, NULL, 'm'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};


int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    //TreeReader tr;
    char * outf = NULL;
    char * treef = NULL;
    char * mrcaf = NULL;
    bool outfileset = false;
    bool fileset = false;
    bool mrcaset = false;
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:o:m:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'm':
                mrcaset = true;
                mrcaf = strdup(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
            case 'V':
                std::cout << versionline << std::endl;
                exit(0);
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(treef, outf);
    }
    
    if (!mrcaset) {
        std::cerr << "Because no file was provided, all the internal nodes" << std::endl;
        std::cerr << "will be labeled" << std::endl;
    }
    
    std::istream * pios = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ofstream * ofstr = NULL;
    
    if (outfileset == true) {
        ofstr = new std::ofstream(outf, std::ios::app);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    if (fileset == true) {
        fstr = new std::ifstream(treef);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }
    
    /* 
       collect clade names
       expecting (new) format:
       MRCANAME = tip1 tip2 ... 
    */
    std::map<std::string, std::vector<std::string> > mrcas;
    if (mrcaset) {
        std::ifstream inmrca(mrcaf);
        std::string mrcaline;
        while (getline(inmrca, mrcaline)) {
            if (mrcaline.empty()) {
                continue;
            }
            std::vector<std::string> searchtokens;
            tokenize(mrcaline, searchtokens, "=");
            std::string mrcaname = searchtokens[0];
            trim_spaces(mrcaname);
            searchtokens.erase(searchtokens.begin());
            searchtokens = tokenize(searchtokens[0]);
            mrcas[mrcaname] = searchtokens;
        }
        inmrca.close();
    }
    
// collect tree(s)
    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        std::cerr << "This really only works with nexus or newick" << std::endl;
        exit(0);
    }
    
    int count = 0;
    bool going = true;
    if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (going) {
                if(mrcaset){
                    std::map<std::string, std::vector<std::string> >::iterator it;
                    for (it = mrcas.begin(); it != mrcas.end(); it++) {
                        //std::cout << "Dealing with clade '" << (*it).first << "'" << std::endl;
                        if (!check_names_against_tree(tree, (*it).second)) {
                            std::cout << "Check mrca file for typos." << std::endl;
                            exit (0);
                        }
                        Node * nd = tree->getMRCA((*it).second);
                        nd->setName((*it).first);
                    }
                } else {
                    for (int i=0; i<tree->getInternalNodeCount();i++) {
                        if (tree->getInternalNode(i)->getName().size() == 0) {
                            tree->getInternalNode(i)->setName("px"+std::to_string(count));
                            count ++;
                        }
                    }
                }
                (*poos) << getNewickString(tree) << std::endl;
                delete tree;
            }
        }
    } else {
        std::map<std::string, std::string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != NULL) {
                if (mrcaset) {
                    std::map<std::string, std::vector<std::string> >::iterator it;
                    for (it = mrcas.begin(); it != mrcas.end(); it++) {
                        //std::cout << "Dealing with clade '" << (*it).first << "'" << std::endl;
                        if (!check_names_against_tree(tree, (*it).second)) {
                            std::cout << "Check mrca file for typos." << std::endl;
                            exit (0);
                        }
                        Node * nd = tree->getMRCA((*it).second);
                        nd->setName((*it).first);
                    }
                } else {
                    for (int i=0; i<tree->getInternalNodeCount();i++) {
                        if (tree->getInternalNode(i)->getName().size() == 0) {
                            tree->getInternalNode(i)->setName("px"+std::to_string(count));
                            count ++;
                        }
                    }
                }
                (*poos) << getNewickString(tree) << std::endl;
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
