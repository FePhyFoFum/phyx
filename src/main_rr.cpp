#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <algorithm>
#include <set>
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
    std::cout << "Reroot (or unroot) a tree file and produce a newick." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN." << std::endl;
    std::cout << "Output is written in newick format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxrr [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE     input tree file, STDIN otherwise" << std::endl;
    std::cout << " -g, --outgroups=CSL  outgroup sep by commas (NO SPACES!)" << std::endl;
    std::cout << " -f, --namesf=FILE    outgroups in a file (each on a line)" << std::endl;
    std::cout << " -r, --ranked         turn on ordering of outgroups. will root on first one present" << std::endl;
    std::cout << " -u, --unroot         unroot the tree" << std::endl;
    std::cout << " -s, --silent         do not error if outgroup(s) not found" << std::endl;
    std::cout << " -o, --outf=FILE      output tree file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help           display this help and exit" << std::endl;
    std::cout << " -V, --version        display version and exit" << std::endl;
    std::cout << " -C, --citation       display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxrr 1.3\n";
    vl += "Copyright (C) 2014-2021 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Stephen A. Smith (blackrim), Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"outgroups", required_argument, nullptr, 'g'},
    {"namesf", required_argument, nullptr, 'f'},
    {"ranked", no_argument, nullptr, 'r'},
    {"unroot", no_argument, nullptr, 'u'},
    {"outf", required_argument, nullptr, 'o'},
    {"silent", no_argument, nullptr, 's'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outgroupsset = false;
    bool namefileset = false;
    bool outfileset = false;
    bool silent = false;
    bool unroot = false;
    bool ranked = false;
    std::vector<std::string> outgroups;

    char * treef = nullptr;
    char * outf = nullptr;
    char * outgroupsc = nullptr;
    char * namesfc = nullptr;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:g:f:ruo:shVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'g':
                outgroupsset = true;
                outgroupsc = strdup(optarg);
                break;
            case 'f':
                namefileset = true;
                namesfc = strdup(optarg);
                check_file_exists(namesfc);
                break;
            case 'r':
                ranked = true;
                break;
            case 'u':
                unroot = true;
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
    
    if (outgroupsset) {
        std::vector<std::string> tokens2;
        tokenize(outgroupsc, tokens2, ",");
        for (auto & tk : tokens2) {
            trim_spaces(tk); // this will never have to be used, as spaces would break cmd line call
            outgroups.push_back(tk);
        }
    } else if (namefileset) {
        std::ifstream nfstr(namesfc);
        std::string tline;
        while (getline_safe(nfstr, tline)) {
            trim_spaces(tline);
            if (!tline.empty()) {
                outgroups.push_back(tline);
            }
        }
        nfstr.close();
        if (!outgroups.empty()) {
            outgroupsset = true;
        } else {
            // empty file
            std::cerr << "Error: no names found. Exiting." << std::endl;
            exit(0);
        }
    }
    if (!outgroupsset && !unroot) {
        std::cerr << "Error: you need to set the outgroup (-g). Exiting." << std::endl;
        exit(0);
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
    
    //read trees 
    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        std::cerr << "Error: this really only works with nexus or newick. Exiting." << std::endl;
        exit(0);
    }
    bool going = true;
    if (!unroot) {
        if (ft == 0) { // Nexus
            std::map<std::string, std::string> translation_table;
            bool ttexists;
            ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
            while (going) {
                Tree * tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                    &translation_table, &going);
                if (tree != nullptr) {
                    if (has_root_edge(tree)) {
                        std::cerr << "Error: tree has a root edge, so rerooting is not possible" << std::endl;
                        std::cerr << "The root edge can be removed with the -r option of pxcltr" << std::endl;
                        std::cerr << "Exiting." << std::endl;
                        exit(0);
                    }
                    bool exists = false;
                    if (ranked) {
                        // find first outgroup present in tree
                        bool ogexists = false;
                        for (const auto & name : outgroups) {
                            if (check_name_against_tree(tree, name)) {
                                std::vector<std::string> og;
                                og.push_back(name);
                                exists = reroot(tree, og, silent);
                                ogexists = true;
                                break;
                            }
                        }
                        // if no valid outgroups, let silent option figure out
                        if (!ogexists) {
                            exists = reroot(tree, outgroups, silent);
                        }
                    } else {
                        exists = reroot(tree, outgroups, silent);
                    }
                    if (!exists) {
                        std::cerr << "The outgroup taxa don't exist in this tree." << std::endl;
                    } else {
                        (*poos) << getNewickString(tree) << std::endl;
                    }
                    delete tree;
                }
            }
        } else if (ft == 1) { // newick
            while (going) {
                Tree * tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
                if (going) {
                    if (has_root_edge(tree)) {
                        std::cerr << "Error: tree has a root edge, so rerooting is not possible" << std::endl;
                        std::cerr << "The root edge can be removed with the -r option of pxcltr" << std::endl;
                        std::cerr << "Exiting." << std::endl;
                        exit(0);
                    }
                    bool exists = false;
                    if (ranked) {
                        // find first outgroup present in tree
                        bool ogexists = false;
                        for (const auto & name : outgroups) {
                            if (check_name_against_tree(tree, name)) {
                                std::vector<std::string> og;
                                og.push_back(name);
                                exists = reroot(tree, og, silent);
                                ogexists = true;
                                break;
                            }
                        }
                        // if no valid outgroups, let silent option figure out
                        if (!ogexists) {
                            exists = reroot(tree, outgroups, silent);
                        }
                    } else {
                        exists = reroot(tree, outgroups, silent);
                    }
                    if (!exists) {
                        std::cerr << "The outgroup taxa don't exist in this tree." << std::endl;
                    } else {
                        (*poos) << getNewickString(tree) << std::endl;
                    }
                    delete tree;
                }
            }
        }
    } else {
        // unroot trees
        if (ft == 0) {
            std::map<std::string, std::string> translation_table;
            bool ttexists;
            ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
            while (going) {
                Tree * tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                    &translation_table, &going);
                if (tree != nullptr) {
                    if (has_root_edge(tree)) {
                        std::cerr << "Error: tree has a root edge, so rerooting is not possible" << std::endl;
                        std::cerr << "The root edge can be removed with the -r option of pxcltr" << std::endl;
                        std::cerr << "Exiting." << std::endl;
                        exit(0);
                    }
                    tree->unRoot();
                    (*poos) << getNewickString(tree) << std::endl;
                    delete tree;
                }
            }
        } else if (ft == 1) {
            while (going) {
                Tree * tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
                if (going) {
                    if (has_root_edge(tree)) {
                        std::cerr << "Error: tree has a root edge, so rerooting is not possible" << std::endl;
                        std::cerr << "The root edge can be removed with the -r option of pxcltr" << std::endl;
                        std::cerr << "Exiting." << std::endl;
                        exit(0);
                    }
                    tree->unRoot();
                    (*poos) << getNewickString(tree) << std::endl;
                    delete tree;
                }
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
