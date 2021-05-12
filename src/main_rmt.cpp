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
#include "constants.h" // contains PHYX_CITATION


void print_help() {
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

std::string versionline("pxrmt 1.2\nCopyright (C) 2014-2021 FePhyFoFum\nLicense GPLv3\nWritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"names", required_argument, NULL, 'n'},
    {"namesf", required_argument, NULL, 'f'},
    {"regex", required_argument, NULL, 'r'},
    {"outf", required_argument, NULL, 'o'},
    {"comp", no_argument, NULL, 'c'},
    {"silent", required_argument, NULL, 's'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
    {NULL, 0, NULL, 0}
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
    std::string regex_pattern = "";

    char * treef = NULL;
    char * outf = NULL;
    char * namesc = NULL;
    char * namesfc = NULL;
    while (1) {
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
    
    if (namesset == true) {
        std::vector<std::string> tokens2;
        std::string del2(",");
        tokens2.clear();
        tokenize(namesc, tokens2, del2);
        for (unsigned int j = 0; j < tokens2.size(); j++) {
            trim_spaces(tokens2[j]);
            names.push_back(tokens2[j]);
        }
    } else if (namefileset == true) {
        std::ifstream nfstr(namesfc);
        std::string tline;
        while (getline_safe(nfstr, tline)) {
            trim_spaces(tline);
            names.push_back(tline);
        }
        nfstr.close();
    } else if (!regex) {
        std::cerr << "Error: you must specify which tips to remove." << std::endl;
        std::cerr << "This can be done with a list (-n) or file (-f) of names, or a regular expression (-r)." << std::endl;
        std::cerr << "Exiting." << std::endl;
        exit(0);
    }

    std::istream * pios = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ofstream * ofstr = NULL;
    
    if (fileset == true) {
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
        poos = &std::cout;
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
    int MAX_RMT = 50;
    
    bool going = true;
    int numLeaves;
    int num_names = 0;
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
            if (going == true) {
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
                
                num_names = currNames.size();
                
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
            if (going == true) {
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
                
                num_names = currNames.size();
                
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
