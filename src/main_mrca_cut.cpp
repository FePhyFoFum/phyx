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
#include "constants.h"

extern std::string PHYX_CITATION;


void print_help() {
    std::cout << "Extract subclade(s) from tree(s)." << std::endl;
    std::cout << "Takes in newick tree and MRCA file with format:" << std::endl;
    std::cout << "MRCANAME = tip1 tip2 ..." << std::endl;
    std::cout << "If multiple MRCAs are provided, multiple subtrees are returned" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxmrcacut [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -t, --treef=FILE    input newick tree file, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output newick file, stout otherwise" << std::endl;
    std::cout << " -m, --mrca=FILE     file containing MRCA declarations" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxmrcacut 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim), Joseph W. Brown");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"mrca", required_argument, NULL, 'm'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    TreeReader tr;
    char * outf = NULL;
    char * treef = NULL;
    char * mrcaf = NULL;
    bool outfileset = false;
    bool fileset = false;
    bool mrcaset = false;
    
    while (1) {
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
    std::ifstream inmrca(mrcaf);
    std::string mrcaline;
    std::map<std::string, std::vector<std::string> > mrcas;
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
    
    // this does not use the conventional tree reader functions
    // update this
    
    // collect tree(s)
    std::vector<std::string> lines;
    std::string line;
    while (getline(*pios, line)) {
        lines.push_back(line);
    }
    
    for (unsigned int i = 0; i < lines.size(); i++) {
        Tree * tree = tr.readTree(lines[i]);
        //std::cout << tree->getExternalNodeCount() << std::endl;

        std::map<std::string, std::vector<std::string> >::iterator it;
        for (it = mrcas.begin(); it != mrcas.end(); it++) {
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
