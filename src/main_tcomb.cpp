#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>
#include <algorithm>
#include <iterator>
#include <cmath>

#include "tree_reader.h"
#include "tree.h"
#include "tree_utils.h"
#include "utils.h"
#include "log.h"
#include "citations.h"


/*
Give two options:
1. pass in 2 trees
2. pass in 1 distribution of trees
*/

void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Combine a set of trees from one file into a tree from another." << std::endl;
    std::cout << "Pass in 2 trees with `t` and `a`." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxtcomb [OPTIONS]... FILE" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE    reference treefile, STDIN otherwise" << std::endl;
    std::cout << " -a, --addtree=FILE  alternate treefile" << std::endl;
    std::cout << " -o, --outf=FILE     output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxtcomb 1.3.2\n";
    vl += "Copyright (C) 2017-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"addtree", required_argument, nullptr, 'a'},
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
    bool addfileset = false;
    
    char * treef = nullptr;
    char * addtreef = nullptr;
    char * outf = nullptr;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:a:o:x:hVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'a':
                addfileset = true;
                addtreef = strdup(optarg);
                check_file_exists(addtreef);
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
    if (addfileset && outfileset) {
        check_inout_streams_identical(addtreef, outf);
    }
    
    std::istream * pios = nullptr;
    std::istream * apios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ifstream * afstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
    if (tfileset) {
        fstr = new std::ifstream(treef);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (!check_for_input_to_stream()) {
            // if both inputs missing: print help
            if (!addfileset) {
                print_help();
                exit(1);
            } else {
                std::cerr << "Error: you need to set an tfile (-t). Exiting." << std::endl;
            exit(1);
            }
        }
    }
    
    if (addfileset) {
        afstr = new std::ifstream(addtreef);
        apios = afstr;
    } else {
        std::cerr << "Error: you need to set an addfile (-a). Exiting." << std::endl;
        exit(0);
    }
    
    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        std::cerr << "Error: this really only works with nexus or newick. Exiting." << std::endl;
        exit(0);
    }
    
    bool going = true;
    Tree * bigtree = nullptr;
    while (going) {
        if (retstring.size() > 1) {
            bigtree = read_next_tree_from_stream_newick (*pios, retstring, &going);
        } else {
            going = false;
        }
    }
    std::set<std::string> btns = bigtree->getRoot()->get_leave_names_set();

    ft = test_tree_filetype_stream(*apios, retstring);
    if (ft != 0 && ft != 1) {
        std::cerr << "Error: this really only works with nexus or newick. Exiting." << std::endl;
        exit(0);
    }
    
    going = true;
    while (going) {
        if (retstring.size() > 1) {
            Tree * addtree = read_next_tree_from_stream_newick (*apios, retstring, &going);
            std::set<std::string> atns = addtree->getRoot()->get_leave_names_set();
            std::vector<std::string> v(btns.size());
            auto it = set_difference(btns.begin(), btns.end(), atns.begin(),
                    atns.end(), v.begin());
            v.resize(static_cast<size_t>(it-v.begin()));
            std::map<std::string, Node *> diffnds;
            std::vector<std::string> diffnms; 
            for (unsigned int i = 0; i < bigtree->getExtantNodeCount(); i++) {
                if (std::count(v.begin(), v.end(), bigtree->getExternalNode(i)->getName()) == 1) {
                    diffnms.push_back(bigtree->getExternalNode(i)->getName());
                }
            } 
            // get start node on big tree by getting the mrca
            Node * connecthere = bigtree->getMRCA(diffnms);
            std::vector<Node *> childs = connecthere->getChildren();
            for (auto & child : childs) {
                std::vector<std::string> v_int(diffnms.size());
                std::set<std::string> lvsset = child->get_leave_names_set();
                auto it1 = set_intersection(lvsset.begin(), lvsset.end(),
                                atns.begin(), atns.end(), v_int.begin());
                v_int.resize(static_cast<size_t>(it1-v_int.begin()));
                if (!v_int.empty()) {
                    // need to add those missing not the overlapping
                    std::vector<std::string> v2(lvsset.size());
                    auto it2 = set_difference(lvsset.begin(), lvsset.end(),
                            atns.begin(), atns.end(), v2.begin());
                    v2.resize(static_cast<size_t>(it2-v2.begin()));
                    for (const auto & v2j : v2) {
                        std::cout << "to add " << v2j << std::endl; 
                        diffnds[v2j] = bigtree->getExternalNode(v2j);
                    }
                    connecthere->removeChild(*child);
                }
            }
            Node * oldroot = addtree->getRoot();
            oldroot->setParent(*connecthere);
            connecthere->addChild(*oldroot);
            addtree->setRoot(connecthere);
            bool didit = false;
            while (!diffnds.empty()) {
                std::cout << "diffnds.size() = " << diffnds.size() << std::endl;
                for (auto iter = diffnds.begin(); iter != diffnds.end(); ++it) {
                    std::cout << iter->first << std::endl;
                    Node * cn = iter->second;
                    while (true) {
                        Node * prn = cn->getParent();
                        std::set<std::string> cs = cn->get_leave_names_set();
                        std::vector<std::string> v_int;
                        set_intersection(cs.begin(), cs.end(), atns.begin(), atns.end(), back_inserter(v_int));
                        if (!v_int.empty()) {
                            std::cout << "this is what we need to add " << prn->getNewick(false) << std::endl;
                            // get nodes
                            // get mrca
                            Node * nd = addtree->getMRCA(v_int);
                            std::cout << "would add to " << nd->getNewick(false) << std::endl;
                            // START HERE
                            std::cout << "v_int.size() = " << v_int.size() << std::endl;
                            if (v_int.size() == 1) {
                                
                                // missing something?
                                
                            }
                            std::vector<std::string> lvsnms = prn->get_leave_names();
                            for (const auto & lvsnm : lvsnms) {
                                if (diffnds.count(lvsnm) > 0) {
                                    diffnds.erase(lvsnm);
                                }
                            }
                            going = false;
                            didit = true;
                            break;
                        }
                        cn = prn;
                    }
                    if (didit) {
                        break;
                    }
                }
            }
            (*poos) << getNewickString(addtree) << std::endl;
            delete addtree;
        } else {
            going = false;
        }
    }

    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
