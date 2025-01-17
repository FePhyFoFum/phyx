
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

#include "utils.h"
#include "tree_reader.h"
#include "relabel.h"
#include "tree_utils.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Taxon relabelling for trees." << std::endl;
    std::cout << "Two ordered lists of taxa, -c (current) and -n (new) must be provided." << std::endl;
    std::cout << "Alternatively, a regex pattern (-p) and replacement (-r) text can be provided." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN." << std::endl;
    std::cout << "Output is written in newick format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxrlt [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE     input tree file, STDIN otherwise" << std::endl;
    std::cout << " -c, --cnames=FILE    file containing current taxon labels (one per line)" << std::endl;
    std::cout << " -n, --nnames=FILE    file containing new taxon labels (one per line)" << std::endl;
    std::cout << " -p, --pattern=STRING regex pattern to replace" << std::endl;
    std::cout << " -r, --replace=STRING replacement pattern" << std::endl;
    std::cout << " -v, --verbose        make the output more verbose" << std::endl;
    std::cout << " -o, --outf=FILE      output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help           display this help and exit" << std::endl;
    std::cout << " -V, --version        display version and exit" << std::endl;
    std::cout << " -C, --citation       display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxrlt 1.3.2\n";
    vl += "Copyright (C) 2018-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph W. Brown, Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"cnames", required_argument, nullptr, 'c'},
    {"nnames", required_argument, nullptr, 'n'},
    {"pattern", required_argument, nullptr, 'p'},
    {"replace", required_argument, nullptr, 'r'},
    {"outf", required_argument, nullptr, 'o'},
    {"verbose", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool tfileset = false;
    bool cfileset = false;
    bool nfileset = false;
    bool verbose = false;
    bool regex = false;
    std::string regex_pattern;
    std::string replacement_text;
    char * outf = nullptr;
    char * treef = nullptr;
    std::string cnamef;
    std::string nnamef;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:c:n:p:r:o:vhVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'c':
                cfileset = true;
                cnamef = strdup(optarg);
                check_file_exists(cnamef);
                break;
            case 'n':
                nfileset = true;
                nnamef = strdup(optarg);
                check_file_exists(nnamef);
                break;
            case 'p':
                regex = true;
                regex_pattern = strdup(optarg);
                break;
            case 'r':
                regex = true;
                replacement_text = strdup(optarg);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'v':
                verbose = true;
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
    
    if (!regex) {
        if (!nfileset || !cfileset) {
            std::cerr << "Error: must supply both name files (-c for current, -n for new). Exiting." << std::endl;
            exit(0);
        }
        
        Relabel rl (cnamef, nnamef, verbose);

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
                    rl.relabel_tree(tree);
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
                    rl.relabel_tree(tree);
                    (*poos) << getNewickString(tree) << std::endl;
                    delete tree;
                }
            }
        }
    } else {
        // regex
        if (replacement_text.empty() || regex_pattern.empty()) {
            std::cerr << "Error: must supply both pattern to match and replacement text. Exiting." << std::endl;
            exit(0);
        }
        
        Relabel rl (regex_pattern, replacement_text);
        
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
                    rl.regex_relabel_tree(tree);
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
                    rl.regex_relabel_tree(tree);
                    (*poos) << getNewickString(tree) << std::endl;
                    delete tree;
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
