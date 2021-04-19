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
#include "constants.h"

extern std::string PHYX_CITATION;


void print_help () {
    std::cout << "Fit a lineage diversification model to a tree." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxbdfit [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE    input treefile, STDIN otherwise" << std::endl;
    std::cout << " -m, --model=STRING  diversification model; either 'yule', 'bd' (default), or 'best'" << std::endl;
    std::cout << " -o, --outf=FILE     output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxbdfit 1.2\nCopyright (C) 2016-2021 FePhyFoFum\nLicense GPLv3\nWritten by Joseph W. Brown");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"model", required_argument, NULL, 'm'},
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
    
    char * treef = NULL;
    char * outf = NULL;
    
    std::string model = "bd";
    std::set<std::string> avail_models{"bd", "yule", "best"};
    
    while (1) {
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
    
    std::istream * pios = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ofstream * ofstr = NULL;

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
            tree = read_next_tree_from_stream_newick (*pios, retstring, &going);
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
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != NULL) {
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
