#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <algorithm>
#include <map>

#include "tree.h"
#include "tree_reader.h"
#include "utils.h"
#include "tree_utils.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Remove tree tips by label." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN." << std::endl;
    std::cout << "Output is written in newick format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxrmt [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE    input tree file, STDIN otherwise" << std::endl;
    std::cout << " -n, --names=CSL     names sep by commas (NO SPACES!)" << std::endl;
    std::cout << " -f, --namesf=FILE   names in a file (each on a line)" << std::endl;
    std::cout << " -r, --regex=STRING  match tip labels by a regular expression" << std::endl;
    std::cout << " -c, --comp          take the complement (i.e. remove any taxa not in list)" << std::endl;
    std::cout << " -s, --silent        suppress warnings of missing tips" << std::endl;
    std::cout << " -o, --outf=FILE     output tree file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxrmt 1.3.2\n";
    vl += "Copyright (C) 2014-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph W. Brown, Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"names", required_argument, nullptr, 'n'},
    {"namesf", required_argument, nullptr, 'f'},
    {"regex", required_argument, nullptr, 'r'},
    {"outf", required_argument, nullptr, 'o'},
    {"comp", no_argument, nullptr, 'c'},
    {"silent", required_argument, nullptr, 's'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
      
    bool fileset = false;
    bool namesset = false;
    bool namefileset = false;
    bool outfileset = false;
    bool silent = false;
    std::vector<std::string> names;
    bool complement = false;
    bool regex = false;
    std::string regex_pattern;

    char * treef = nullptr;
    char * outf = nullptr;
    char * namesc = nullptr;
    char * namesfc = nullptr;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:n:f:r:co:shVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'n':
                namesset = true;
                namesc = strdup(optarg);
                break;
            case 'f':
                namefileset = true;
                namesfc = strdup(optarg);
                check_file_exists(namesfc);
                break;
            case 'r':
                regex = true;
                regex_pattern = strdup(optarg);
                break;
            case 'c':
                complement = true;
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 's':
                silent = true;
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
    
    if (fileset && outfileset) {
        check_inout_streams_identical(treef, outf);
    }
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
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
    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    if (namesset) {
        std::vector<std::string> tokens2;
        std::string del2(",");
        tokens2.clear();
        tokenize(namesc, tokens2, del2);
        for (auto & tk : tokens2) {
            trim_spaces(tk); // this will never have to be used, as spaces would break cmd line call
            names.push_back(tk);
        }
    } else if (namefileset) {
        std::ifstream nfstr(namesfc);
        std::string tline;
        while (getline_safe(nfstr, tline)) {
            trim_spaces(tline);
            names.push_back(tline);
        }
        nfstr.close();
    } else if (!regex) {
        std::cerr << "Error: you must specify which tips to remove." << std::endl;
        std::cerr << "This can be done with a list (-n) or file (-f) of names, or a regular expression (-r)."
                << std::endl;
        std::cerr << "Exiting." << std::endl;
        exit(0);
    }
    
    //read trees 
    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        std::cerr << "Error: this really only works with nexus or newick. Exiting." << std::endl;
        exit(0);
    }
    
    // magic number: if number to be remove is > than this, use the induced tree procedure instead
    // this is demonstrated here: https://github.com/FePhyFoFum/phyx/issues/74
    unsigned int MAX_RMT = 50;
    
    bool going = true;
    unsigned int numLeaves;
    unsigned int num_names = 0;
    std::vector<std::string> toKeep; // if trace is used
    std::vector<std::string> currNames; // keep original copy of names in case sampling differs across trees
    
    if (ft == 0) {
        // nexus
        std::map<std::string, std::string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (going) {
                numLeaves = tree->getExternalNodeCount();
                currNames = names; 
                
                if (regex) {
                    currNames = get_names_in_tree_regex(tree, regex_pattern);
                }
                if (complement) {
                    currNames = get_complement_tip_set(tree, currNames);
                }
                // check names against the tree (instead of just length of names, which might be bad)
                if (!regex && !complement) {
                    currNames = get_names_in_tree(tree, currNames);
                }
                
                num_names = static_cast<unsigned int>(currNames.size());
                
                if (num_names == 0) {
                    if (!silent) {
                        std::cerr << "Error: no matching tip labels. Returning original tree." << std::endl;
                    }
                    (*poos) << getNewickString(tree) << std::endl;
                } else if (numLeaves - num_names < 2) {
                    std::cerr << "Error: pruning would produce a tree with "
                        << (numLeaves - num_names) << " tips. No result is returned." << std::endl;
                } else {
                    if (num_names < MAX_RMT) {
                        remove_tips(tree, currNames, silent);
                    } else {
                        toKeep = get_complement_tip_set(tree, currNames);
                        tree = get_induced_tree(tree, toKeep, silent);
                    }
                    (*poos) << getNewickString(tree) << std::endl;
                }
                delete tree;
            }
        }
    } else if (ft == 1) {
        // newick
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (going) {
                numLeaves = tree->getExternalNodeCount();
                currNames = names; 
                
                if (regex) {
                    currNames = get_names_in_tree_regex(tree, regex_pattern);
                }
                if (complement) {
                    currNames = get_complement_tip_set(tree, currNames);
                }
                // check names against the tree (instead of just length of names, which might be bad)
                if (!regex && !complement) {
                    currNames = get_names_in_tree(tree, currNames);
                }
                
                num_names = static_cast<unsigned int>(currNames.size());
                
                if (num_names == 0) {
                    if (!silent) {
                        std::cerr << "Error: no matching tip labels. Returning original tree." << std::endl;
                    }
                    (*poos) << getNewickString(tree) << std::endl;
                } else if (numLeaves - num_names < 2) {
                    std::cerr << "Error: pruning would produce a tree with "
                        << (numLeaves - num_names) << " tips. No result is returned." << std::endl;
                }  else {
                    if (num_names < MAX_RMT) {
                        remove_tips(tree, currNames, silent);
                    } else {
                        toKeep = get_complement_tip_set(tree, currNames);
                        tree = get_induced_tree(tree, toKeep, silent);
                    }
                    (*poos) << getNewickString(tree) << std::endl;
                }
                delete tree;
            }
        }
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
