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
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
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

std::string get_version_line () {
    std::string vl = "pxcltr 1.3.2\n";
    vl += "Copyright (C) 2017-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"root", no_argument, nullptr, 'r'},
    {"labels", no_argument, nullptr, 'l'},
    {"knuckles", no_argument, nullptr, 'k'},
    {"outf", required_argument, nullptr, 'o'},
    {"showd", no_argument, nullptr, 's'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
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
    
    char * treef = nullptr;
    char * outf = nullptr;
    
    while (true) {
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
                std::cout << get_version_line() << std::endl;
                exit(0);
            case 'C':
                std::cout << get_phyx_citation() << std::endl;
                exit(0);
            default:
                print_error(*argv);
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
    std::istream* pios = nullptr;
    std::ostream* poos = nullptr;
    std::ifstream* fstr = nullptr;
    std::ofstream* ofstr = nullptr;

    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    if (tfileset) {
        fstr = new std::ifstream(treef);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (!check_for_input_to_stream()) {
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
        while (going) {
            Tree * tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (tree != nullptr) {
                CleanTree ct(tree, removerootedge, removelabels, removeknuckles);
                (*poos) << getNewickString(tree) << std::endl;
                delete tree;
            }
        }
    } else if (ft == 0) { // Nexus. need to worry about possible translation tables
        std::map<std::string, std::string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        while (going) {
            Tree * tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != nullptr) {
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
