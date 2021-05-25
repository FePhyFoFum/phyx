#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <algorithm>
#include <set>

#include "tree.h"
#include "tree_reader.h"
#include "utils.h"
#include "log.h"
#include "citations.h"


void print_help (void);
std::string get_version_line (void);

void print_help () {
    std::cout << "This will print out bipartitions found in treefile." << std::endl;
    std::cout << "Trees are assumed rooted unless the -e argument is provided." << std::endl;
    std::cout << "This will take a newick- or nexus-formatted tree from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxbp [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE    input treefile, STDIN otherwise" << std::endl;
    std::cout << " -v, --verbose       give more output" << std::endl;
    std::cout << " -e, --edgeall       force edgewise (not node - so when things are unrooted) and" << std::endl; 
    std::cout << "                           assume all taxa are present in all trees" << std::endl;
    std::cout << " -u, --uniquetree    output unique trees and *no* other output" << std::endl;
    std::cout << " -m, --maptree=FILE  put the bipart freq on the edges of this tree. This will " << std::endl;
    std::cout << "                           create a *.pxbpmapped.tre file." << std::endl;
    std::cout << " -c, --cutoff        skip biparts that have support lower than this." << std::endl;
    std::cout << " -s, --suppress      don't print all the output (maybe you use this" << std::endl;
    std::cout << "                           with the maptree feature" << std::endl;
    std::cout << " -o, --outf=FILE     output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxbp 1.2\n";
    vl += "Copyright (C) 2017-2021 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"verbose", no_argument, nullptr, 'v'},
    {"edgeall", no_argument, nullptr, 'e'},
    {"uniquetree", no_argument, nullptr, 'u'},
    {"maptree", required_argument, nullptr, 'm'},
    {"cutoff", required_argument, nullptr, 'c'},
    {"suppress", no_argument, nullptr, 's'},
    {"first", no_argument, nullptr, 'f'},
    {"outf", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    bool mapfileset = false;
    bool verbose = false;
    bool edgewisealltaxa = false;
    bool uniquetree = false;
    bool suppress = false;
    bool cutoff = false;
    bool firsttree = false;
    char * treef = nullptr;
    char * mtreef = nullptr;
    char * outf = nullptr;
    double cutnum = 0;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:o:m:c:vseufhVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'v':
                verbose = true;
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'e':
                edgewisealltaxa = true;
                break;
            case 'u':
                uniquetree = true;
                break;
            case 's':
                suppress = true;
                break;
            case 'f':
                firsttree = true;
                break;
            case 'c':
                cutnum = string_to_double(optarg, "-c");
                cutoff = true;
                break;
            case 'm':
                mapfileset = true;
                mtreef = strdup(optarg);
                check_file_exists(mtreef);
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
    
    if (fileset && outfileset) {
        check_inout_streams_identical(treef, outf);
    }
    
    std::istream * pios = nullptr;
    std::istream * mpios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ifstream * mfstr = nullptr;
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

    // map file
    if (mapfileset) {
        mfstr = new std::ifstream(mtreef);
        mpios = mfstr;
    }

    //------READ TREES
    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        std::cerr << "Error: this really only works with nexus or newick. Exiting." << std::endl;
        exit(0);
    }

    std::vector<Tree *> trees;
    bool going = true;
    if (ft == 0) {
        std::map<std::string, std::string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (tree != nullptr) {
                trees.push_back(tree);
            }
        }
    } else if (ft == 1) {
        Tree * tree;
        while (going) {
            tree = read_next_tree_from_stream_newick(*pios, retstring, &going);
            if (going) {
                trees.push_back(tree);
            } 
        }
    }
    //-----END READ TREES
    //-----READ MAP TREE
    Tree * maptree = nullptr;
    if (mapfileset) {
        ft = test_tree_filetype_stream(*mpios, retstring);
        if (ft != 0 && ft != 1) {
            std::cerr << "Error: this really only works with nexus or newick. Exiting." << std::endl;
            exit(0);
        }

        going = true;
        if (ft == 0) {
            std::map<std::string, std::string> translation_table;
            bool ttexists;
            ttexists = get_nexus_translation_table(*mpios, &translation_table, &retstring);
            while (going) {
                maptree = read_next_tree_from_stream_nexus(*mpios, retstring, ttexists,
                    &translation_table, &going);
                if (maptree != nullptr) {
                    going = false;
                    break;
                }
            }
        } else if (ft == 1) {
            while (going) {
                maptree = read_next_tree_from_stream_newick(*mpios, retstring, &going);
                if (going) {
                    going = false;
                    break;
                } 
            }
        } 
    }
    //----END READ MAP TREE

    int numtrees = static_cast<int>(trees.size());
    if (numtrees == 0) {
        if (fileset) {
            fstr->close();
            delete pios;
        }
        if (outfileset) {
            ofstr->close();
            delete poos;
        }
        std::cerr << "there are no trees;" << std::endl;
    }
    
    // get the biparts for the trees
    std::vector<std::string> names;
    std::set<std::string> names_s;
    std::map<std::string, int> name_index;
    std::map<int, std::string> name_st_index;
    // added to make sure we have all the names if it is partially overlapping
    for (auto & tree : trees) {
        for (int i = 0; i < tree->getExternalNodeCount(); i++) {
            std::string tname = tree->getExternalNode(i)->getName();
            if (find(names.begin(), names.end(), tname) == names.end()) {
                name_index[tname] = i;
                names.push_back(tname);
                name_st_index[i] = tname;
            }
        }
    }
    copy(names.begin(), names.end(), inserter(names_s, names_s.begin()));

    std::vector<std::vector<int> > biparts; // first part of the bipart
    std::vector<std::vector<int> > biparts2; // second part of the bipart
    std::vector<std::vector<int> > not_included; // the names that aren't in the tree
    std::vector<double> bp_count;
    for (unsigned long i = 0; i < static_cast<unsigned long>(numtrees); i++) {
        // get the biparts
        bool unrooted = false;
        int numch = trees[i]->getRoot()->getChildCount();
        if (numch > 2) {
            unrooted = true;
        }
        std::vector<std::string> rt_nms = trees[i]->getRoot()->get_leave_names();
        std::set<std::string> rt_nms_set;
        copy(rt_nms.begin(), rt_nms.end(), inserter(rt_nms_set, rt_nms_set.begin()));
        // get the ones that aren't in the tree at all
        std::vector<std::string> not_included_nms(rt_nms_set.size());
        std::vector<int> not_included_i;
        std::vector<std::string>::iterator it2;
        it2 = std::set_difference(names_s.begin(), names_s.end(), rt_nms_set.begin(),
                rt_nms_set.end(), not_included_nms.begin());
        not_included_nms.resize(it2-not_included_nms.begin());
        for (const auto & not_included_nm : not_included_nms) {
            std::cerr << " not included: " << not_included_nm << std::endl;
            not_included_i.push_back(name_index[not_included_nm]);
        }
        std::vector<int> bp_count_tree; // for edgewise to make sure we don't double count
        for (int j = 0; j < trees[i]->getInternalNodeCount(); j++) {
            std::vector<std::string> nms = trees[i]->getInternalNode(j)->get_leave_names();
            // skip the root
            if (nms.size() == rt_nms.size()) {
                continue;
            }
            // if we are using a cutoff, skip the edge that is below the num
            if (cutoff) {
                if (trees[i]->getInternalNode(j)->getName().length() < 1) {
                    continue;
                }
                char* pEnd;
                double td = strtod(trees[i]->getInternalNode(j)->getName().c_str(),&pEnd);
                if (td < cutnum) {
                    continue;
                }
            }
            // end using cutoffs
            std::vector<int> nms_i;
            std::set<std::string> nms_s;
            copy(nms.begin(), nms.end(), inserter(nms_s, nms_s.begin()));
            for (const auto & nm : nms) {
                nms_i.push_back(name_index[nm]);
            }
            sort(nms_i.begin(), nms_i.end());
            // get the other side of the bipart
            std::vector<int> nms_i2;
            std::vector<std::string> nms_s2(rt_nms.size());
            std::vector<std::string>::iterator it;
            it = std::set_difference(rt_nms_set.begin(), rt_nms_set.end(), nms_s.begin(),
                nms_s.end(), nms_s2.begin());
            nms_s2.resize(it-nms_s2.begin());
            for (const auto & k : nms_s2) {
                nms_i2.push_back(name_index[k]);
            }
            sort(nms_i2.begin(), nms_i2.end());
            // check to see if the bipart is new
            if (!edgewisealltaxa) {
                // this is nodewise and we dont' assume we have all the taxa
                if (static_cast<int>(std::count(biparts.begin(), biparts.end(), nms_i)) == 0 && 
                    static_cast<int>(std::count(biparts2.begin(), biparts2.end(), nms_i2)) == 0) {
                    biparts.push_back(nms_i);
                    biparts2.push_back(nms_i2);
                    not_included.push_back(not_included_i);
                    bp_count.push_back(1);
                } else {
                    // get index 
                    size_t index = find(biparts.begin(), biparts.end(), nms_i)-biparts.begin();
                    bp_count[index] += 1;
                }
                // do the otherside for unrooted
                if (unrooted && trees[i]->getInternalNode(j)->getParent() == trees[i]->getRoot()) {
                    if (static_cast<int>(std::count(biparts.begin(), biparts.end(), nms_i2)) == 0 && 
                        static_cast<int>(std::count(biparts2.begin(), biparts2.end(), nms_i)) == 0) {
                        biparts.push_back(nms_i2);
                        biparts2.push_back(nms_i);
                        not_included.push_back(not_included_i);
                        bp_count.push_back(1);
                    } else {
                        size_t index = find(biparts.begin(), biparts.end(), nms_i2)-biparts.begin();
                        bp_count[index] += 1;
                    }
                }
            } else {
                // this is edgewise and we assume all the taxa 
                // this is for reporting a b | c d instead of a b and c d separately
                
                // first check to make sure that both sides at least have two taxa
                if (nms_i.size() < 2 || nms_i2.size() < 2) {
                    continue;
                }
                if (static_cast<int>(std::count(biparts.begin(), biparts.end(), nms_i)) == 0 && 
                    static_cast<int>(std::count(biparts.begin(), biparts.end(), nms_i2)) == 0) {
                    biparts.push_back(nms_i);
                    biparts2.push_back(nms_i2);
                    not_included.push_back(not_included_i);
                    bp_count.push_back(1);
                    bp_count_tree.push_back(bp_count.size()-1);
                } else {
                    // get index 
                    size_t index;
                    if (static_cast<int>(std::count(biparts.begin(), biparts.end(), nms_i)) == 1) {
                        index = find(biparts.begin(), biparts.end(), nms_i)-biparts.begin();
                    } else {
                        index = find(biparts.begin(), biparts.end(), nms_i2)-biparts.begin();
                    }
                    // need to accommodate that it can be reflected in the edgewise case
                    if (static_cast<int>(std::count(bp_count_tree.begin(), bp_count_tree.end(), index)) == 0) {
                        bp_count[index] += 1;
                        bp_count_tree.push_back(index);
                    }
                }

            }
        }
    }
    if (!suppress) {
        (*poos) << numtrees << " trees " << std::endl;
        (*poos) << biparts.size() << " unique clades found" << std::endl;
        // calculate the logical matrix of biparts for each tree
        // the matrix will have each i as a tree and 
        //   each matrix[i] will represent a 1 if the tree has the bipart and 0 if not
        int colsize = static_cast<int>(biparts.size()) + 1;
        if (edgewisealltaxa) {
            colsize--; // no root for edgewise
        }
        std::vector<int> cols(static_cast<unsigned long>(colsize), 0);
        std::vector<std::vector<int> > matrix(static_cast<unsigned long>(numtrees), cols);
        for (unsigned long i = 0; i < static_cast<unsigned long>(numtrees); i++) {
            bool unrooted = false;
            int numch = trees[i]->getRoot()->getChildCount();
            if (numch > 2) {
                unrooted = true;
            }
            for (int j = 0; j < trees[i]->getInternalNodeCount(); j++) {
                // if we are using a cutoff, skip the edge that is below the num
                if (cutoff) {
                    if (trees[i]->getInternalNode(j)->getName().length() < 1) {
                        continue;
                    }
                    char* pEnd;
                    double td = strtod(trees[i]->getInternalNode(j)->getName().c_str(),&pEnd);
                    if (td < cutnum) {
                        continue;
                    }
                }
                // end using cutoffs
                // should this do the root? I don't think so
                if (edgewisealltaxa && trees[i]->getInternalNode(j) == trees[i]->getRoot()) {
                    continue;
                }

                std::vector<std::string> nms = trees[i]->getInternalNode(j)->get_leave_names();
                std::vector<int> nms_i;
                for (const auto & nm : nms) {
                    nms_i.push_back(name_index[nm]);
                }
                sort(nms_i.begin(), nms_i.end());
                size_t index;
                if (edgewisealltaxa) {
                    if (static_cast<int>(std::count(biparts.begin(), biparts.end(), nms_i)) == 1) {
                        index = find(biparts.begin(), biparts.end(), nms_i)-biparts.begin();
                    } else {
                        index = find(biparts2.begin(), biparts2.end(), nms_i)-biparts2.begin();
                    } 
                } else {
                    index = find(biparts.begin(), biparts.end(), nms_i) - biparts.begin();
                }
                matrix[i][index] = 1;
                if (unrooted && trees[i]->getInternalNode(j)->getParent() == trees[i]->getRoot()) {
                    std::vector<std::string> rt_nms = trees[i]->getRoot()->get_leave_names();
                    std::set<std::string> rt_nms_set;
                    copy(rt_nms.begin(), rt_nms.end(), inserter(rt_nms_set, rt_nms_set.begin()));
                    std::set<std::string> nms_s;
                    copy(nms.begin(), nms.end(), inserter(nms_s, nms_s.begin()));
                    nms_i.clear();
                    std::vector<int> nms_i2;
                    std::vector<std::string> nms_s2(rt_nms.size());
                    std::vector<std::string>::iterator it;
                    it = std::set_difference(rt_nms_set.begin(), rt_nms_set.end(), nms_s.begin(),
                            nms_s.end(), nms_s2.begin());
                    nms_s2.resize(it-nms_s2.begin());
                    for (const auto & k : nms_s2) {
                        nms_i2.push_back(name_index[k]);
                    }
                    for (const auto & ii : nms_s) {
                        nms_i.push_back(name_index[ii]);
                    }
                    sort(nms_i2.begin(), nms_i2.end());
                    sort(nms_i.begin(), nms_i.end());
                    int x = find(biparts.begin(), biparts.end(), nms_i2) - biparts.begin();
                    if (x == static_cast<int>(matrix[i].size())) {
                        x = find(biparts.begin(), biparts.end(), nms_i)-biparts.begin();
                    }
                    matrix[i][static_cast<unsigned long>(x)] = 1;
                }
            }
        }
        // print the unique trees
        if (uniquetree) {
            std::cout << "====UNIQUE TREES====" << std::endl;
            std::set<std::string> un_trees;
            for (unsigned long i = 0; i < static_cast<unsigned long>(numtrees); i++) {
                std::string sv = get_string_vector(matrix[i]);
                if (count(un_trees.begin(), un_trees.end(), sv) == 0) {
                    std::cout << trees[i]->getRoot()->getNewick(false) << ";" << std::endl;
                    un_trees.insert(sv);
                }
            }
            std::cout << "==END UNIQUE TREES==" << std::endl;
            exit(0);
        }

        // constructing the logical matrix
        // the logical matrix has each row as a bipart1 and each col as a name
        // there is a one if the bipart has the name
        // need to add the -1 for bipart2
        std::vector<int> cols2(names.size(), 0);
        std::vector<std::vector<int> > logical_matrix (biparts.size(), cols2);
        for (unsigned int i = 0; i < biparts.size(); i++) {
            for (unsigned int j = 0; j < names.size(); j++) {
                if (std::count(biparts[i].begin(), biparts[i].end(), name_index[names[j]]) != 0) {
                    logical_matrix[i][j] = 1;
                }
            }
            //std::cout << get_string_vector(logical_matrix[i]) << std::endl;
        }
        
        double smallest_proportion = 0.0;
        double TSCA = 0;
        // get the conflicting bipartitions
        // initialize results vectors
        for (unsigned int i = 0; i < biparts.size(); i++) {
            if (firsttree) {
                if (matrix[0][i] != 1) {
                    continue;
                }
            }
            unsigned int sumc = sum_matrix_col(matrix, i);
            if (sumc != trees.size() && sumc > (smallest_proportion*trees.size())) {
                std::vector<std::string> nms;
                for (unsigned int k = 0; k < biparts[i].size(); k++) {
                    nms.push_back(name_st_index[biparts[i][k]]);
                }
                (*poos) << "CLADE: " << get_string_vector(nms);
                if (edgewisealltaxa) {
                    std::vector<std::string> nms_o;
                    for (unsigned int k = 0; k < biparts2[i].size(); k++) {
                        nms_o.push_back(name_st_index[biparts2[i][k]]);
                    }
                    (*poos) << "| " << get_string_vector(nms_o);
                }
                double totalcount = bp_count[i];
                std::vector<double> conflict_nums;
                conflict_nums.push_back(bp_count[i]);
                if (verbose) {
                    (*poos) << "\n\tCONFLICTS:" << std::endl;
                }
                for (unsigned int j = 0; j < biparts.size(); j++) {
                    unsigned int sumc2 = sum_matrix_col(matrix, j);
                    if (i != j && sumc2 != trees.size() && sumc2 > (smallest_proportion*trees.size())) {
                        bool logitest = test_logical(logical_matrix[i], logical_matrix[j], edgewisealltaxa);
                        if (logitest) {
                            std::vector<std::string> nms2;
                            for (unsigned int k = 0; k < biparts[j].size(); k++) {
                                nms2.push_back(name_st_index[biparts[j][k]]);
                            }    
                            totalcount += bp_count[j];
                            conflict_nums.push_back(bp_count[j]);
                            if (verbose) {
                                (*poos) << " \t " << get_string_vector(nms2);
                                if (edgewisealltaxa) {
                                    std::vector<std::string> nms_o;
                                    for (unsigned int k = 0; k < biparts2[j].size(); k++) {
                                        nms_o.push_back(name_st_index[biparts2[j][k]]);
                                    }
                                    (*poos) << "| " << get_string_vector(nms_o);
                                }
                                (*poos) << "\tCOUNT:\t" << bp_count[j] << "\tTREEFREQ:\t"
                                        << bp_count[j]/trees.size() << std::endl;
                            }
                        }
                    }
                }
                // calculate ICA
                double sign = 1;
                for (unsigned int j = 0; j < conflict_nums.size(); j++) {
                    conflict_nums[j] /= totalcount;
                    if (conflict_nums[j] > conflict_nums[0]) {
                        sign = -1;
                    }
                }
                double ICA = 1; // same as logn(conflict_nums.size(), conflict_nums.size());
                for (unsigned int j = 0; j < conflict_nums.size(); j++) {
                    ICA += (conflict_nums[j] * logn(conflict_nums[j], conflict_nums.size()));
                }
                TSCA += ICA;
                ICA *= sign;
                (*poos) << "\tFREQ:\t" << conflict_nums[0] << "\tICA:\t" << ICA << "\tCOUNT:\t"
                        << bp_count[i] << "\tTREEFREQ:\t" << bp_count[i]/trees.size() << std::endl;
                if (verbose) {
                    (*poos) << "\tTREES:\t";
                    for (unsigned int j = 0; j < matrix.size(); j++) {
                        if (matrix[j][i] == 1) {
                            (*poos) << j << " ";   
                        }
                    }
                    (*poos) << std::endl;
                }
            } else if (sumc == trees.size()) {
                std::vector<std::string> nms;
                for (unsigned int k = 0; k < biparts[i].size(); k++) {
                    nms.push_back(name_st_index[biparts[i][k]]);
                }
                (*poos) << "CLADE: " << get_string_vector(nms);
                if (edgewisealltaxa) {
                    std::vector<std::string> nms_o;
                    for (unsigned int k = 0; k < biparts2[i].size(); k++) {
                        nms_o.push_back(name_st_index[biparts2[i][k]]);
                    }
                    (*poos) << "| " << get_string_vector(nms_o);
                }
                (*poos) << "\tFREQ:\t1.\tICA:\t1.\tCOUNT:\t" << bp_count[i] << "\tTREEFREQ:\t1." << std::endl;
                TSCA += 1;
                if (verbose) {
                    (*poos) << "\tTREES:\t";
                    for (unsigned int j = 0; j < matrix.size(); j++) {
                        if (matrix[j][i] == 1) {
                            (*poos) << j << " ";   
                        }
                    }
                    (*poos) << std::endl;
                }
            }
        }
        (*poos) << "TSCA: " << TSCA << std::endl;
    } // end suppress if
    if (mapfileset) {
        std::string mot(mtreef);
        mot = mot +".pxbpmapped.tre";
        std::ofstream * mofstr = new std::ofstream(mot);
        std::ostream * mpoos = mofstr;
        for (int i = 0; i < maptree->getInternalNodeCount(); i++) {
            if (maptree->getInternalNode(i) == maptree->getRoot()) {
                continue;
            }
            std::vector<std::string> nms = maptree->getInternalNode(i)->get_leave_names();
            std::vector<int> nms_i;
            for (const auto & nm : nms) {
                nms_i.push_back(name_index[nm]);
            }
            sort(nms_i.begin(), nms_i.end());
            size_t index;
            bool found = false;
            if (edgewisealltaxa) {
                if (static_cast<int>(std::count(biparts.begin(), biparts.end(), nms_i)) == 1) {
                    index = find(biparts.begin(), biparts.end(), nms_i)-biparts.begin();
                    found = true;
                } else {
                    index = find(biparts2.begin(), biparts2.end(), nms_i)-biparts2.begin();
                    found = true;
                } 
            } else {
                index = find(biparts.begin(), biparts.end(), nms_i) - biparts.begin();
                found = true;
            }
            if (!found) {
                maptree->getInternalNode(i)->setName("0.0");
            } else {
                maptree->getInternalNode(i)->setName(std::to_string(bp_count[index]/trees.size()));
            }
        }
        (*mpoos) << maptree->getRoot()->getNewick(true) << ";" << std::endl;
        mofstr->close();
        delete mpoos;
    }

    // shut things down
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
