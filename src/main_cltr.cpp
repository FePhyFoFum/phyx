#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>
#include <algorithm>
#include <cmath>

#include "tree_reader.h"
#include "tree.h"
#include "tree_utils.h"
#include "utils.h"
#include "clean_tree.h"
#include "log.h"
#include "constants.h" // contains PHYX_CITATION


void print_help () {
    std::cout << "General tree cleaner." << std::endl;
    std::cout << "Removes annotations (node labels), 'knuckles' (2-degree nodes), and" << std::endl;
    std::cout << "root edges to generate a 'vanilla' newick representation." << std::endl;
    std::cout << "By default removes all properties. Alternatively choose 1 property." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN." << std::endl;
    std::cout << "Output is written in newick format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxcltr [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE    input treefile, STDIN otherwise" << std::endl;
    std::cout << " -r, --root          remove root edge (if present)" << std::endl;
    std::cout << " -l, --labels        remove internal node labels" << std::endl;
    std::cout << " -k, --knuckles      remove 2-degree internal nodes" << std::endl;
    std::cout << " -o, --outf=FILE     output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxcltr 1.2\nCopyright (C) 2017-2021 FePhyFoFum\nLicense GPLv3\nWritten by Joseph W. Brown");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"root", no_argument, NULL, 'r'},
    {"labels", no_argument, NULL, 'l'},
    {"knuckles", no_argument, NULL, 'k'},
    {"outf", required_argument, NULL, 'o'},
    {"showd", no_argument, NULL, 's'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool tfileset = false;
    bool optionsset = false;
    bool removerootedge = false;
    bool removelabels = false;
    bool removeknuckles = false;
    
    // need option to write nexus
    
    char * treef = NULL;
    char * outf = NULL;
    
    while(true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:rlko:x:hVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'r':
                removerootedge = true;
                optionsset = true;
                break;
            case 'l':
                removelabels = true;
                optionsset = true;
                break;
             case 'k':
                removeknuckles = true;
                optionsset = true;
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
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
    
    if (tfileset && outfileset) {
        check_inout_streams_identical(treef, outf);
    }
    
    // if no options set, do all by default
    if (!optionsset) {
        removerootedge = true;
        removelabels = true;
        removeknuckles = true;    }
    /*
    if ((removerootedge + removelabels) > 1) {
        std::cerr << "Error: specify 1 property only (or leave blank to clean all). Exiting." << std::endl;
        exit(0);
    }
    */
    std::istream* pios = NULL;
    std::ostream* poos = NULL;
    std::ifstream* fstr = NULL;
    std::ofstream* ofstr = NULL;

    if (outfileset == true) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    if (tfileset == true) {
        fstr = new std::ifstream(treef);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }
    
    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        std::cerr << "Error: this really only works with nexus or newick. Exiting." << std::endl;
        exit(0);
    }
    
    bool going = true;
    if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (tree != NULL) {
                CleanTree ct(tree, removerootedge, removelabels, removeknuckles);
                (*poos) << getNewickString(tree) << std::endl;
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
            if (tree != NULL) {
                CleanTree ct(tree, removerootedge, removelabels, removeknuckles);
                (*poos) << getNewickString(tree) << std::endl;
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
