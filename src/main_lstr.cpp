#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstring>
#include <getopt.h>

#include "tree_info.h"
#include "tree_reader.h"
#include "tree.h"
#include "utils.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Print tree summary." << std::endl;
    std::cout << "By default returns all properties. Alternatively choose 1 property." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxlstr [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE    input tree file, STDIN otherwise" << std::endl;
    std::cout << " -r, --rooted        return whether the tree is rooted" << std::endl;
    std::cout << " -a, --age           return the height of root (must be rooted and ultrametric)" << std::endl;
    std::cout << " -n, --ntips         return the number of terminals" << std::endl;
    std::cout << " -u, --ultrametric   return whether tree is ultrametric" << std::endl;
    std::cout << " -b, --binary        return whether tree is binary" << std::endl;
    std::cout << " -l, --length        return the length of the tree" << std::endl;
    std::cout << " -i, --tiplabels     return all tip labels (one per line)" << std::endl;
    std::cout << " -v, --rtvar         return root-to-tip variance" << std::endl;
    std::cout << " -o, --outf=FILE     output tree stats file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxlstr 1.3.1\n";
    vl += "Copyright (C) 2016-2024 FePhyFoFum\n";
    vl += "License: GPLv3\n";
    vl += "Written by Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"rooted", no_argument, nullptr, 'r'},
    {"age", no_argument, nullptr, 'a'},
    {"ntips", no_argument, nullptr, 'n'},
    {"ultrametric", no_argument, nullptr, 'u'},
    {"binary", no_argument, nullptr, 'b'},
    {"length", no_argument, nullptr, 'l'},
    {"tiplabels", no_argument, nullptr, 'i'},
    {"rtvar", no_argument, nullptr, 'v'},
    {"outf", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    bool optionsset = false; // if true, do not return all properties
    int propcount = 0; // count how many properties are requested
    bool ultracheck = false;
    bool binarycheck = false;
    bool lengthcheck = false;
    bool agecheck = false;
    bool rootedcheck = false;
    bool ntipcheck = false;
    bool namecheck = false;
    bool rtvarcheck = false;
    char * outf = nullptr;
    char * treef = nullptr;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:vranublio:x:hVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'r':
                rootedcheck = true;
                optionsset = true;
                propcount++;
                break;
            case 'a':
                agecheck = true;
                optionsset = true;
                propcount++;
                break;
            case 'n':
                ntipcheck = true;
                optionsset = true;
                propcount++;
                break;
            case 'u':
                ultracheck = true;
                optionsset = true;
                propcount++;
                break;
            case 'b':
                binarycheck = true;
                optionsset = true;
                propcount++;
                break;
            case 'l':
                lengthcheck = true;
                optionsset = true;
                propcount++;
                break;
            case 'i':
                namecheck = true;
                optionsset = true;
                propcount++;
                break;
            case 'v':
                rtvarcheck = true;
                optionsset = true;
                propcount++;
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
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
    if (outfileset) {
        ofstr = new std::ofstream(outf);
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
    
    if (fileset && outfileset) {
        check_inout_streams_identical(treef, outf);
    }
    
    if (propcount > 1) {
        std::cerr << "Error: specify 1 property only (or leave blank to show all properties). Exiting." << std::endl;
        exit(0);
    }

    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        std::cerr << "Error: this really only works with nexus or newick. Exiting." << std::endl;
        exit(0);
    }
    
    int treeCounter = 0;
    bool going = true;
    if (ft == 1) {
        while (going) {
            Tree * tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (tree != nullptr) {
                if (!optionsset) {
                    (*poos) << "tree #: " << treeCounter << std::endl;
                    TreeInfo ti(tree);
                    ti.get_stats(poos);
                    delete tree;
                    treeCounter++;
                } else {
                    // only a single property
                    TreeInfo ti(tree, ultracheck, binarycheck, agecheck, rootedcheck,
                        ntipcheck, lengthcheck, namecheck, rtvarcheck, poos);
                }
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
                if (!optionsset) {
                    (*poos) << "tree #: " << treeCounter << std::endl;
                    TreeInfo ti(tree);
                    ti.get_stats(poos);
                    delete tree;
                    treeCounter++;
                } else {
                    // only a single property
                    TreeInfo ti(tree, ultracheck, binarycheck, agecheck, rootedcheck,
                        ntipcheck, lengthcheck, namecheck, rtvarcheck, poos);
                }
            }
        }
    }
    
    return EXIT_SUCCESS;
}
