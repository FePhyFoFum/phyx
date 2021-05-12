#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <map>

#include "tree.h"
#include "tree_reader.h"
#include "utils.h"
#include "tree_utils.h"
#include "log.h"
#include "citations.h" // contains PHYX_CITATION


void print_help() {
    std::cout << "This will convert a tree file to vanilla Nexus format." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxt2nex [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE    input tree file, STDIN otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output tree file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxt2nex 1.2\nCopyright (C) 2021 FePhyFoFum\nLicense GPLv3\nWritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"outf", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    char * treef = nullptr;
    char * outf = nullptr;
    while(true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:o:hVC", long_options, &oi);
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
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
    if (fileset ) {
        fstr = new std::ifstream(treef);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (!check_for_input_to_stream()) {
            print_help();
            exit(1);
        }
    }
    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    //read trees 
    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    
    int treeCounter = 0;
    bool going = true;
    
    if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (tree != nullptr) {
                if (treeCounter == 0) {
                    (*poos) << "#NEXUS" << std::endl << "Begin trees;" << std::endl;
                }
                if (is_rooted(tree)) {
                    (*poos) << "tree tree" << treeCounter << " = [&R] "
                            << getNewickString(tree) << std::endl;
                } else {
                    (*poos) << "tree tree" << treeCounter << " = [&U] "
                            << getNewickString(tree) << std::endl;
                }
                treeCounter++;
            }
        }
        (*poos) << "end;" << std::endl;
    } else if (ft == 0) { // Nexus. need to worry about possible translation tables
        std::map<std::string, std::string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != nullptr) {
                if (treeCounter == 0) {
                    (*poos) << "#NEXUS" << std::endl << "Begin trees;" << std::endl;
                }
                if (is_rooted(tree)) {
                    (*poos) << "tree tree" << treeCounter << " = [&R] "
                            << getNewickString(tree) << std::endl;
                } else {
                    (*poos) << "tree tree" << treeCounter << " = [&U] "
                            << getNewickString(tree) << std::endl;
                }
                treeCounter++;
            }
        }
        (*poos) << "end;" << std::endl;
    } else {
        std::cerr << "Error: tree format not recognized. Exiting." << std::endl;
        exit(1);
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
