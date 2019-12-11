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
#include "constants.h"

extern std::string PHYX_CITATION;


/*
Give two options:
1. pass in 2 trees
2. pass in 1 distribution of trees
*/

void print_help () {
    std::cout << "Combine a set of trees from one file into a tree from another." << std::endl;
    std::cout << "Pass in 2 trees with `t` and `a`." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxtcomb [OPTIONS]... FILE" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE    reference treefile, stdin otherwise" << std::endl;
    std::cout << " -a, --addtree=FILE  alternate treefile" << std::endl;
    std::cout << " -o, --outf=FILE     output file, stout otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxtcomb 1.0\nCopyright (C) 2017-2020 FePhyFoFum\nLicense GPLv3\nWritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"treef", required_argument, NULL, 't'},
    {"addtree", required_argument, NULL, 'a'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool tfileset = false;
    bool addfileset = false;
    
    char * treef = NULL;
    char * addtreef = NULL;
    char * outf = NULL;
    
    while (1) {
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
    if (addfileset && outfileset) {
        check_inout_streams_identical(addtreef, outf);
    }
    
    std::istream * pios = NULL;
    std::istream * apios = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ifstream * afstr = NULL;
    std::ofstream * ofstr = NULL;

    if (outfileset == true) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    if (addfileset == true) {
        afstr = new std::ifstream(addtreef);
        apios = afstr;
    } else {
        std::cerr << "Error: you need to set an addfile (-a). Exiting." << std::endl;
        exit(0);
    }
    
    if (tfileset == true) {
        fstr = new std::ifstream(treef);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (check_for_input_to_stream() == false) {
            std::cerr << "Error: you need to set an tfile (-t). Exiting." << std::endl;
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
    Tree * bigtree = NULL;
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
            std::vector<std::string>::iterator it = set_difference(btns.begin(), btns.end(), atns.begin(), atns.end(), v.begin());
            v.resize(it-v.begin());
            std::map<std::string, Node *> diffnds;
            std::vector<std::string> diffnms; 
            for (int i = 0; i <bigtree->getExtantNodeCount(); i++) {
                if (count(v.begin(), v.end(), bigtree->getExternalNode(i)->getName())==1) {
                    diffnms.push_back(bigtree->getExternalNode(i)->getName());
                }
            } 
            //get start node on big tree by getting the mrca
            Node * connecthere = bigtree->getMRCA(diffnms);
            std::vector<Node *> childs = connecthere->getChildren();
            for (unsigned int i=0; i < childs.size(); i++) {
                std::vector<std::string> v_int(diffnms.size());
                std::vector<std::string>::iterator it;
                std::set<std::string> lvsset = childs[i]->get_leave_names_set();
                it = set_intersection(lvsset.begin(), lvsset.end(),
                                atns.begin(), atns.end(), v_int.begin());
                v_int.resize(it-v_int.begin());
                if (v_int.size() > 0) {
                    //need to add those missing not the overlapping
                    std::vector<std::string> v2(lvsset.size());
                    std::vector<std::string>::iterator it2 = set_difference(lvsset.begin(), lvsset.end(), atns.begin(), atns.end(), v2.begin());
                    v2.resize(it2-v2.begin());
                    for(unsigned int j=0;j<v2.size();j++) {
                        std::cout << "to add " <<  v2[j]<< std::endl; 
                        diffnds[v2[j]] = bigtree->getExternalNode(v2[j]);
                    }
                    connecthere->removeChild(*childs[i]);
                }
            }
            Node * oldroot = addtree->getRoot();
            oldroot->setParent(*connecthere);
            connecthere->addChild(*oldroot);
            addtree->setRoot(connecthere);
            bool didit = false;
            while (diffnds.size() > 0) {
                for (std::map<std::string, Node*>::iterator it = diffnds.begin(); it != diffnds.end(); ++it) {
                    //std::cout << it->first << std::endl;
                    Node * cn = it->second;
                    bool goi = true;
                    while (goi == true) {
                        Node * prn = cn->getParent();
                        std::set<std::string> cs = cn->get_leave_names_set();
                        std::vector<std::string> v_int;
                        set_intersection(cs.begin(), cs.end(), atns.begin(), atns.end(), back_inserter(v_int));
                        if (v_int.size() > 0) {
                            std::cout <<  "this is what we need to add " << prn->getNewick(false) << std::endl;
                            //get nodes
                            //get mrca
                            Node * nd = addtree->getMRCA(v_int);
                            std::cout << "would add to " <<  nd->getNewick(false) << std::endl;
                            //START HERE
                            if (v_int.size() == 1) {
                                
                            }
                            std::vector<std::string> lvsnms = prn->get_leave_names();
                            for(unsigned int i=0; i < lvsnms.size(); i++) {
                                if (diffnds.count(lvsnms[i]) > 0) {
                                    diffnds.erase(lvsnms[i]);
                                }
                            }
                            going = false;
                            didit = true;
                            break;
                        }
                        cn = prn;
                    }
                    if (didit == true) {
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
