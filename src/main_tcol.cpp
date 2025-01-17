#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <map>
#include <getopt.h>

#include "utils.h"
#include "tree_reader.h"
#include "tree_utils.h"
#include "log.h"
#include "aa2cdn.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Add information to a tree so that you can color the edges." << std::endl;
    std::cout << "This will take nexus and newick inputs from a file or STDIN." << std::endl;
    std::cout << "Results are written in nexus format so that it can be read by figtree." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxtcol [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE     input tree file, STDIN otherwise" << std::endl;
    std::cout << " -m, --mrcaf=FILE     file with mrcas and annotations, tab separated" << std::endl;
    std::cout << " -d, --nodeidf=FILE   file with nodeids (labels) and annotations, tab separated" << std::endl;
    std::cout << " -o, --outf=FILE      output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help           display this help and exit" << std::endl;
    std::cout << " -V, --version        display version and exit" << std::endl;
    std::cout << " -C, --citation       display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxtcol 1.3.2\n";
    vl += "Copyright (C) 2016-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"mrcaf", required_argument, nullptr, 's'},
    {"nodeidf", required_argument, nullptr, 'r'},
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
    
    // hrm add this back in?
    //bool mrcafset = false; not used
    //char * mrcaf = nullptr;
    
    bool nodeidfset = false;
    char * nodeidf = nullptr;
    char * outf = nullptr;
    char * treef = nullptr;

    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:m:d:o:hVC", long_options, &oi);
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
                //mrcafset = true;
                //mrcaf = strdup(optarg);
                break;
            case 'd':
                nodeidfset = true;
                nodeidf = strdup(optarg);
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
    
    //for node ids
    std::map<std::string, std::string> nodeid_map;
    if (nodeidfset) {
        std::ifstream nfstr(nodeidf);
        std::string tline;
        while (getline_safe(nfstr, tline)) {
            trim_spaces(tline);
            if (tline.empty()) {
                continue;
            }
            std::vector<std::string> tokens2;
            tokenize(tline, tokens2, "\t");
            for (auto & tk : tokens2) {
                trim_spaces(tk);
            }
            if (tokens2.size() != 2) {
                continue;
            }
            nodeid_map[tokens2[0]] = tokens2[1];
        }
        nfstr.close();
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
        (*poos) << "#NEXUS" << std::endl;
        (*poos) << "begin trees;" << std::endl;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (going) {
                for (unsigned int i = 0; i < tree->getInternalNodeCount(); i++) {
                    Node * tnode = tree->getInternalNode(i);
                    if (nodeid_map.find(tnode->getName()) != nodeid_map.end()) {
                        tnode->setName("[&name=\""+tnode->getName()+"\",ann="+nodeid_map[tnode->getName()]+"]");
                    } 
                }
                for (unsigned int i = 0; i < tree->getExternalNodeCount(); i++) {
                    Node * tnode = tree->getExternalNode(i);
                    if (nodeid_map.find(tnode->getName()) != nodeid_map.end()) {
                        tnode->setName(tnode->getName()+"[&ann="+nodeid_map[tnode->getName()]+"]");
                    } 
                }
                // put annotations here
                (*poos) << "tree tree = " << getNewickString(tree) << std::endl;
                delete tree;
            }
        }
        (*poos) << "end;" << std::endl;
    } else if (ft == 0) { // Nexus. need to worry about possible translation tables
        std::map<std::string, std::string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        Tree * tree;
        (*poos) << "#NEXUS" << std::endl;
        (*poos) << "begin trees;" << std::endl;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != nullptr) {
                for (unsigned int i = 0; i < tree->getInternalNodeCount(); i++) {
                    Node * tnode = tree->getInternalNode(i);
                    if (nodeid_map.find(tnode->getName()) != nodeid_map.end()) {
                        tnode->setName("[&name=\""+tnode->getName()+"\",ann="+nodeid_map[tnode->getName()]+"]");
                    } 
                }
                for (unsigned int i = 0; i < tree->getExternalNodeCount(); i++) {
                    Node * tnode = tree->getExternalNode(i);
                    if (nodeid_map.find(tnode->getName()) != nodeid_map.end()) {
                        tnode->setName(tnode->getName()+"[&ann="+nodeid_map[tnode->getName()]+"]");
                    } 
                }   
                // put annotations here
                (*poos) << "tree tree = " << getNewickString(tree) << std::endl;
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
