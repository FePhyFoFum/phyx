#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "tree_utils.h"
#include "log.h"
#include "citations.h" // contains PHYX_CITATION


void print_help() {
    std::cout << "Extract subclade(s) from tree(s)." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN," << std::endl;
    std::cout << "and an MRCA file with format:" << std::endl;
    std::cout << "MRCANAME = tip1 tip2 ..." << std::endl;
    std::cout << "If multiple MRCAs are provided, multiple subtrees are returned" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxmrcacut [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE    input newick tree file, STDIN otherwise" << std::endl;
    std::cout << " -m, --mrca=FILE     file containing MRCA declarations" << std::endl;
    std::cout << " -o, --outf=FILE     output newick file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxmrcacut 1.2\nCopyright (C) 2013-2021 FePhyFoFum\nLicense GPLv3\nWritten by Stephen A. Smith (blackrim), Joseph W. Brown");

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"outf", required_argument, nullptr, 'o'},
    {"mrca", required_argument, nullptr, 'm'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    TreeReader tr;
    char * outf = nullptr;
    char * treef = nullptr;
    char * mrcaf = nullptr;
    bool outfileset = false;
    bool fileset = false;
    bool mrcaset = false;
    
    while(true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:o:m:hVC", long_options, &oi);
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
            case 'C':
                std::cout << PHYX_CITATION << std::endl;
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
        std::cerr << "Error: must supply mrca file. Exiting." << std::endl;
        exit(0);
    }
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
    if (outfileset) {
        ofstr = new std::ofstream(outf, std::ios::app);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    if (fileset) {
        fstr = new std::ifstream(treef);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (!check_for_input_to_stream()) {
            print_help();
            exit(1);
        }
    }
    
    /* 
       collect clade names
       expecting (new) format:
       MRCANAME = tip1 tip2 ... 
    */
    std::ifstream inmrca(mrcaf);
    std::string mrcaline;
    std::map<std::string, std::vector<std::string> > mrcas;
    while (getline_safe(inmrca, mrcaline)) {
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
    
    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        std::cerr << "Error: this really only works with nexus or newick. Exiting." << std::endl;
        exit(0);
    }
    
    bool going = true;
    std::map<std::string, std::vector<std::string> >::iterator it;
    
    if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (tree != nullptr) {
                for (it=mrcas.begin(); it != mrcas.end(); it++) {
                    //std::cout << "Dealing with clade '" << (*it).first << "'" << std::endl;
                    if (!check_names_against_tree(tree, (*it).second)) {
                        // allow more flexibility here
                        std::cerr << "Error: check mrca file for typos. Exiting." << std::endl;
                        exit(0);
                    }
                    Node * nd = tree->getMRCA((*it).second);
                    bool bl = has_branchlengths(tree);
                    (*poos) << nd->getNewick(bl) << ";" << std::endl;
                }
                delete tree;
            }
        }
    } else if (ft == 0) { // Nexus. need to worry about possible translation tables
        std::map<std::string, std::string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != nullptr) {
                for (it=mrcas.begin(); it != mrcas.end(); it++) {
                    //std::cout << "Dealing with clade '" << (*it).first << "'" << std::endl;
                    if (!check_names_against_tree(tree, (*it).second)) {
                        // allow more flexibility here
                        std::cerr << "Error: check mrca file for typos. Exiting." << std::endl;
                        exit(0);
                    }
                    Node * nd = tree->getMRCA((*it).second);
                    bool bl = has_branchlengths(tree);
                    (*poos) << nd->getNewick(bl) << ";" << std::endl;
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
