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
#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "log.h"

void print_help() {
    std::cout << "This will print out partitions found in seqfile." << std::endl;
    std::cout << "Currently only supports DNA alignments. << std::endl";
    std::cout << "Can read from stdin or file." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxbpsq [OPTION]... [FILE]..." << std::endl;
    std::cout << std::endl;
    std::cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << std::endl;
    std::cout << " -t, --treef=FILE    input tree file, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output file, stout otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}
/*
 * add you name if you contribute (probably add another line)
 */
std::string versionline("pxbpsq 0.1\nCopyright (C) 2014 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim), Joseph W. Brown");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool tfileset = false;
    bool outfileset = false;
    bool conditional = false; // not used
    char * seqf = NULL;
    char * treef = NULL;
    char * outf = NULL;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:t:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 't':
                tfileset = true;
                conditional = true;
                treef = strdup(optarg);
                check_file_exists(treef);
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
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
    }
    
    std::istream * pios = NULL;
    std::istream * piost = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ofstream * ofstr = NULL;
    std::ifstream * tfstr = NULL;
    
    if (fileset == true) {
        fstr = new std::ifstream(seqf);
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
    // if treefile there, read the trees
    // then we will calculate the conditional biparts as well
    std::vector<Tree *> trees;
    if (tfileset == true) {
        tfstr = new std::ifstream(treef);
        piost = tfstr;
        TreeReader tr;
        std::string retstring;
        while (getline(*piost,retstring)) {
            trees.push_back(tr.readTree(retstring));
        }
        int numtrees = trees.size();    
        if (numtrees == 0) {
            std::cout << "there are no trees;" << std::endl;
        }
    }

    // read the sequences
    std::vector<Sequence> seqs;
    Sequence seq;
    std::string retstring;
    
    int ft = test_seq_filetype_stream(*pios,retstring);
    while (read_next_seq_from_stream(*pios,ft,retstring,seq)) {
        seqs.push_back(seq);
    }
    // fasta has a trailing one
    if (ft == 2) {
        seqs.push_back(seq);
    }
    
    // get the biparts for the trees
    std::vector<std::string> names;
    std::map<std::string, int> name_index;
    std::map<int, std::string> name_st_index;
    for (unsigned int i=0; i < seqs.size(); i++) {
        std::string tname = seqs[i].get_id();
        name_index[tname] = i;
        (*poos) << tname << " " << i << std::endl;
        names.push_back(tname);
        name_st_index[i] = tname;
    }
    int numseqs = seqs.size();
    int numcols = seqs[0].get_sequence().size();

    // hardcoded for DNA
    // can try to get more elegant later
    std::vector<std::vector<std::vector<int> > > all_bp; 
    std::vector<double> bp_count;
    std::vector<int> seq_bipart_map;
    // for each column
    for (int i=0; i < numcols; i++) {
        std::vector<int> tbpa;
        std::vector<int> tbpc;
        std::vector<int> tbpg;
        std::vector<int> tbpt;
        int skip = 0;
        for (int j=0; j < numseqs; j++) {
            if (seqs[j].get_sequence()[i] == 'A') {
                tbpa.push_back(j);
            } else if (seqs[j].get_sequence()[i] == 'C') {
                tbpc.push_back(j);
            } else if (seqs[j].get_sequence()[i] == 'G') {
                tbpg.push_back(j);
            } else if (seqs[j].get_sequence()[i] == 'T') {
                tbpt.push_back(j);
            } else {
                std::cout << "don't recognize character " << seqs[j].get_sequence()[i] << std::endl;
                skip += 1;
            }
        }
        if (skip == numseqs) {
            continue;
        }
        sort(tbpa.begin(),tbpa.end());
        sort(tbpc.begin(),tbpc.end());
        sort(tbpg.begin(),tbpg.end());
        sort(tbpt.begin(),tbpt.end());
        std::vector<std::vector<int> > tallbp;
        tallbp.push_back(tbpa);
        tallbp.push_back(tbpc);
        tallbp.push_back(tbpg);
        tallbp.push_back(tbpt);
        // check to see if the bipart is new
        bool add = true;
        for (unsigned int j=0; j < all_bp.size(); j++) { // for each bipart set
            bool got = true;
            for (unsigned int k=0; k < tallbp.size(); k++) { // or each nucleotide set
                //bool tm = true; not used
                if ((int)count(all_bp[j].begin(),all_bp[j].end(),tallbp[k]) == 0) { // no match
                    //tm = false;
                    got = false;
                    break;
                }
            }
            if (got == true) { // we have a match
                bp_count[j] += 1;
                seq_bipart_map.push_back(j);
                add = false;
                break;
            }
        }
        if (add == true) {
            all_bp.push_back(tallbp);
            bp_count.push_back(1);
            seq_bipart_map.push_back(bp_count.size()-1);
        }
    }
    
    (*poos) << numseqs << " sequences with " << numcols << " cols" <<   std::endl;
    (*poos) << all_bp.size() << " unique parts found" << std::endl;
    
    // calculate the ICA
    // get the parts that are not compatible
    // this would be an intersection of each of the parts and if one is null then it is compatible
    double ACA = 0;
    std::vector<double> bp_ica(bp_count.size());
    for (unsigned int i=0; i < all_bp.size(); i++) {
        double totalcount = bp_count[i];
        std::vector<double> conflict_nums;
        conflict_nums.push_back(bp_count[i]);
        std::vector<int> conflicts;
        for (unsigned int j=0; j < all_bp.size(); j++) {
            if (i == j) {
                continue;
            }
            bool good = true;
            int compcount = 0;
            for (unsigned int m=0; m < all_bp[i].size(); m++) {
                if (all_bp[i][m].size() <= 1) {
                    continue;
                } else {
                    compcount += 1;
                }
                int badcount = 0;
                int compcount2 = 0;
                for (unsigned int n=0; n < all_bp[i].size(); n++) {
                    if (all_bp[j][n].size() <= 1) {
                        continue;
                    }
                    compcount2 += 1;
                    // intersection of [i][m] and [j][n]
                    std::vector<int> v3;
                    set_intersection(all_bp[i][m].begin(),all_bp[i][m].end(),all_bp[j][n].begin(),all_bp[j][n].end(),back_inserter(v3));
                    if (v3.size() > 0 && v3.size() < all_bp[j][n].size() && v3.size() < all_bp[i][m].size()) {
                        badcount += 1;
                    }
                }
                if (badcount >= 2 && compcount > 1 && compcount2 > 1) {
                    good = false;
                    break;
                }
            }
            if (good == false) {
                conflicts.push_back(j);
                conflict_nums.push_back(bp_count[j]);
                totalcount += bp_count[j];
    //        std::cout << i << " ("<<bp_count[i] << ") and " << j << " ("  <<bp_count[j] << " )" << std::endl;
    //        std::cout << "\t" << get_string_vector(all_bp[i][0]) << " | " << get_string_vector(all_bp[i][1]) << " | " << get_string_vector(all_bp[i][2]) << " | " << get_string_vector(all_bp[i][3]) << std::endl;
    //        std::cout << "\t" << get_string_vector(all_bp[j][0]) << " | " << get_string_vector(all_bp[j][1]) << " | " << get_string_vector(all_bp[j][2]) << " | " << get_string_vector(all_bp[j][3]) << std::endl;
            }
        }
        if (conflict_nums.size() == 1) {
            (*poos) << i  << " ("<<bp_count[i] << ") no conflict " << get_string_vector(all_bp[i][0]) << " | " <<  get_string_vector(all_bp[i][1]) << " | " << get_string_vector(all_bp[i][2]) << " | " << get_string_vector(all_bp[i][3]) << std::endl;
            bp_ica[i] = 1;
            continue;
        }
        // calculate ICA
        double sign = 1;
        for (unsigned int j=0; j < conflict_nums.size(); j++) {
            conflict_nums[j]/=totalcount;
            if (conflict_nums[j] > conflict_nums[0]) {
                sign = -1;
            }
        }
        double ICA = 1; // same as logn(conflict_nums.size(),conflict_nums.size());
        for (unsigned int j=0; j < conflict_nums.size(); j++) {
            ICA += (conflict_nums[j]*logn(conflict_nums[j],conflict_nums.size()));
        }
        ACA += ICA;
        ICA *= sign;
        bp_ica[i] = ICA;
        (*poos) << i <<" ("  << bp_count[i] << ")\t" << ICA << " " << get_string_vector(all_bp[i][0]) << " | " <<  get_string_vector(all_bp[i][1]) << " | " << get_string_vector(all_bp[i][2]) << " | " << get_string_vector(all_bp[i][3]) << std::endl;
    }
    (*poos) << " " << ACA << std::endl;
 
    //std::cout << get_string_vector(seq_bipart_map) << std::endl;

    // tree processing
    if (tfileset) {
        for (unsigned int t=0; t < trees.size(); t++) {
            std::vector<std::string> rt_nms = trees[t]->getRoot()->get_leave_names();
            std::set<std::string> rt_nms_set;
            copy(rt_nms.begin(),rt_nms.end(),inserter(rt_nms_set,rt_nms_set.begin()));
            for (int j=0; j < trees[t]->getInternalNodeCount(); j++) {
                std::vector<std::string> nms = trees[t]->getInternalNode(j)->get_leave_names();
                (*poos) << get_string_vector(nms) << std::endl;
                std::vector<int> nms_i;
                std::set<std::string> nms_s;
                copy(nms.begin(),nms.end(),inserter(nms_s,nms_s.begin()));
                for (unsigned int k=0; k < nms.size(); k++) {
                    nms_i.push_back(name_index[nms[k]]);
                }
                sort(nms_i.begin(),nms_i.end());
                // get the other side of the bipart
                std::vector<int> nms_i2;
                std::vector<std::string> nms_s2(rt_nms.size());
                std::vector<std::string>::iterator it;
                it = set_difference(rt_nms_set.begin(),rt_nms_set.end(),nms_s.begin(),nms_s.end(),nms_s2.begin());
                nms_s2.resize(it-nms_s2.begin());
                for (unsigned int k=0; k < nms_s2.size(); k++) {
                    nms_i2.push_back(name_index[nms_s2[k]]);
                }
                // find the biparts that are the same
                std::vector<int> matches;
                for (unsigned int i=0; i < all_bp.size(); i++) {
                    for (unsigned int m=0; m < all_bp[i].size(); m++) {
                        if (all_bp[i][m].size() <= 1 || (all_bp[i][m].size() != nms_i.size() && all_bp[i][m].size() != nms_i2.size() ) ) {
                            continue;
                        }
                        std::vector<int> v3;
                        std::set_intersection(all_bp[i][m].begin(),all_bp[i][m].end(),nms_i.begin(),nms_i.end(),back_inserter(v3));
                        if (v3.size() == nms_i.size()) {
                            (*poos) << i << " (" << bp_ica[i] << ") ";
                            break;
                        }
                        std::vector<int> v4;
                        std::set_intersection(all_bp[i][m].begin(),all_bp[i][m].end(),nms_i2.begin(),nms_i2.end(),back_inserter(v4));
                        if (v4.size() == nms_i2.size()) {
                            (*poos) << i << " (" << bp_ica[i] << ") ";
                            break;
                        }
                    } //each part of the bipart
                } //each bipart
                (*poos) << std::endl;
            } //each internal node
        } //each tree
        tfstr->close();
        delete piost;
    }

    //shut things down
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
