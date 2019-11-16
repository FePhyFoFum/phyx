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

void print_help() {
    std::cout << "This will convert a tree file to vanilla Nexus format." << std::endl;
    std::cout << "Can read from stdin or file." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxt2nex [OPTION]... [FILE]..." << std::endl;
    std::cout << std::endl;
    std::cout << " -t, --treef=FILE    input tree file, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output tree file, stout otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxt2nex 0.1\nCopyright (C) 2019 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    char * treef = NULL;
    char * outf = NULL;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:o:hV", long_options, &oi);
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
            default:
                print_error(argv[0], (char)c);
                exit(0);

        }
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(treef, outf);
    }
    
    std::istream * pios = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ofstream * ofstr = NULL;
    
    if (fileset == true ) {
        fstr = new std::ifstream(treef);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }
    if (outfileset == true) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = & std::cout;
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
            if (tree != NULL) {
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
            if (tree != NULL) {
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
