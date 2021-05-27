#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

#include "utils.h"
#include "tree_reader.h"
#include "tscale.h"
#include "tree_utils.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << "Tree rescaling by providing either scaling factor or root height." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN." << std::endl;
    std::cout << "Output is written in newick format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxtscale [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE        input tree file, STDIN otherwise" << std::endl;
    std::cout << " -s, --scale=DOUBLE      edge length scaling factor" << std::endl;
    std::cout << " -r, --rootheight=DOUBLE height of root (tree must be ultrametric)" << std::endl;
    std::cout << " -o, --outf=FILE         output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help              display this help and exit" << std::endl;
    std::cout << " -V, --version           display version and exit" << std::endl;
    std::cout << " -C, --citation          display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxtscale 1.2\n";
    vl += "Copyright (C) 2016-2021 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"scale", required_argument, nullptr, 's'},
    {"rootheight", required_argument, nullptr, 'r'},
    {"outf", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool tfileset = false;
    double rootheight = 0.0;
    double scalef = 1.0;
    bool heightset = false;
    bool scaleset = false;
    char * outf = nullptr;
    char * treef = nullptr;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:s:r:o:hVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 's':
                scalef = string_to_double(optarg, "-s");
                scaleset = true;
                break;
            case 'r':
                rootheight = string_to_double(optarg, "-r");
                heightset = true;
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
                print_error(argv[0]);
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
    
    if (heightset && scaleset) {
        std::cerr << "Error: supply only rootheight (-r) or scale (-s), not both. Exiting." << std::endl;
        exit(0);
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
    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    TScale ts;
    
    if (heightset) {
        ts.set_rootheight(rootheight);
    } else {
        ts.set_scalef(scalef);
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
            if (going) {
                if (heightset) {
                    // have to check ultrametricity
                    bool isultra = is_ultrametric_paths(tree);
                    if (!isultra) {
                        std::cerr << "Error: setting root height only works for ultrametric trees. Exiting."
                                << std::endl;
                        exit(0);
                    }
                }
                ts.rescale(tree);
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
                if (heightset) {
                    // have to check ultrametricity
                    bool isultra = is_ultrametric_paths(tree);
                    if (!isultra) {
                        std::cerr << "Error: setting root height only works for ultrametric trees. Exiting."
                                << std::endl;
                        exit(0);
                    }
                }
                ts.rescale(tree);
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
