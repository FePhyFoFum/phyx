#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>
#include <algorithm>
#include <cmath>
#include <set>

#include "tree_reader.h"
#include "tree.h"
#include "tree_utils.h"
#include "utils.h"
#include "bd_fit.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Fit a lineage diversification model to a tree." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxbdfit [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE    input treefile, STDIN otherwise" << std::endl;
    std::cout << " -m, --model=STRING  diversification model. either 'yule', 'bd' (default), or 'best'" << std::endl;
    std::cout << " -o, --outf=FILE     output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxbdfit 1.3.2\n";
    vl += "Copyright (C) 2016-2024 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"model", required_argument, nullptr, 'm'},
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
    
    char * treef = nullptr;
    char * outf = nullptr;
    
    std::string model = "bd";
    std::set<std::string> avail_models{"bd", "yule", "best"};
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:m:o:x:hVC", long_options, &oi);
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
                // check valid models here
                model = string_to_lower(optarg);
                if (avail_models.find(model) == avail_models.end()) {
                    std::cerr << "Error: model '" << model << "' is not recognized. Exiting."
                            << std::endl;
                    exit(0);
                }
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
            Tree * tree = read_next_tree_from_stream_newick (*pios, retstring, &going);
            if (going) {
                // in addition to checking ultramtericity, the following sets node heights
                if (is_ultrametric_paths(tree)) {
                    BDFit bd(tree, model);
                    bd.get_pars(poos);
                    delete tree;
                } else {
                    std::cerr << "Error: tree is not ultrametric. Exiting." << std::endl;
                    exit(0);
                }
            }
        }
    }
    if (ft == 0) {
        std::map<std::string, std::string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        while (going) {
            Tree * tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != nullptr) {
                // in addition to checking ultramtericity, the following sets node heights
                if (is_ultrametric_paths(tree)) {
                    BDFit bd(tree, model);
                    bd.get_pars(poos);
                    delete tree;
                } else {
                    std::cerr << "Error: tree is not ultrametric. Exiting." << std::endl;
                    exit(0);
                }
            }
        }
    }
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
