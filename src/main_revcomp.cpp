#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <time.h>
#include <cstdlib>
#include <getopt.h>

#include "seq_reader.h"
#include "sequence.h"
#include "seq_utils.h"
#include "utils.h"
#include "log.h"
#include "edlib.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << "Reverse complement sequences." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << "Results are written in fasta format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxrevcomp [OPTIONS]... [FILE]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -i, --ids=IDS       a comma sep list of ids to flip (NO SPACES!)" << std::endl;
    std::cout << " -g, --guess         EXPERIMENTAL: guess whether there are seqs that need to be " << std::endl;
    std::cout << "                       rev comp. uses edlib library on first seq" << std::endl;
    std::cout << " -p, --pguess        EXPERIMENTAL: progressively guess " << std::endl;
    std::cout << " -m, --sguess        EXPERIMENTAL: sampled guess " << std::endl;
    std::cout << " -o, --outf=FILE     output sequence file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxrevcomp 1.2\n";
    vl += "Copyright (C) 2017-2021 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"ids", required_argument, nullptr, 'i'},
    {"guess", no_argument, nullptr, 'g'},
    {"pguess", no_argument, nullptr, 'p'},
    {"sguess", no_argument, nullptr, 'm'},
    {"outf", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

bool reverse_it_or_not(std::vector<Sequence>& seqs, Sequence comp_seq) {
    int best_distance = 10000000;
    int best_dis_rev = 100000000;
    std::string comp = comp_seq.get_sequence();
    std::string revcomp = comp_seq.reverse_complement();
    for (auto & seq : seqs) {
        EdlibAlignResult result = edlibAlign(comp.c_str(), comp.length(), seq.get_sequence().c_str(), 
                seq.get_length(), edlibNewAlignConfig(best_distance, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
        if (result.editDistance < 0) {
            continue;
        }
        if (result.editDistance < best_distance) {
            best_distance = result.editDistance;
        }
        edlibFreeAlignResult(result);
    }
    for (auto & seq : seqs) {
        EdlibAlignResult result = edlibAlign(revcomp.c_str(), revcomp.length(), seq.get_sequence().c_str(), 
                seq.get_length(), edlibNewAlignConfig(best_distance, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
        if (result.editDistance < 0) {
            continue;
        }
        if (result.editDistance < best_dis_rev) {
            best_dis_rev = result.editDistance;
        }
        edlibFreeAlignResult(result);
    }
    if (best_dis_rev < best_distance) {
        return true;
    }
    return false;
}

int main(int argc, char * argv[]) {
    srand (time(nullptr));
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    bool idsset = false;
    std::vector<std::string> ids;
    
    bool guess = false;
    bool pguess = false;
    bool sguess = false;
    double sguess_samplenum = 0.2; // 10% of them will be used for revcomp
    char * seqf = nullptr;
    char * outf = nullptr;
    char * idssc = nullptr;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:i:o:mgphVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'i':
                idsset = true;
                idssc = strdup(optarg);
                break;
            case 'g':
                guess = true;
                break;
            case 'p':
                guess = true;
                pguess = true;
                break;
            case 'm':
                guess = true;
                sguess = true;
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
                print_error(argv[0]);
                exit(0);
        }
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    if (idsset) {
        std::vector<std::string> tokens2;
        tokenize(idssc, tokens2, ",");
        for (auto & tk : tokens2) {
            trim_spaces(tk); // this will never have to be used, as spaces would break cmd line call
            ids.push_back(tk);
        }
    }

    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
    if (fileset) {
        fstr = new std::ifstream(seqf);
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
    
    Sequence seq;
    std::string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    if (!guess) {
        while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
            if (!idsset || std::count(ids.begin(), ids.end(), seq.get_id()) == 1) {
                seq.perm_reverse_complement();
            }
            (*poos) << seq.get_fasta();
        }
        if (ft == 2) {
            if (!idsset || std::count(ids.begin(), ids.end(), seq.get_id()) == 1) {
                seq.perm_reverse_complement();
            }
            (*poos) << seq.get_fasta();
        }
    } else {
       bool first = true;
       std::vector<Sequence> done; //for pguess
       while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
           if (first) {
               done.push_back(seq);
               (*poos) << seq.get_fasta();
               first = false;
           } else {
               if (reverse_it_or_not(done, seq)) {
                    seq.perm_reverse_complement();
               }
               (*poos) << seq.get_fasta();
               if (pguess) {
                   done.push_back(seq);
               } else if (sguess) {
                   double r = (static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
                    if (r < sguess_samplenum) {
                        done.push_back(seq);
                    }
               }
           }
        }
        if (ft == 2) {
           if (reverse_it_or_not(done, seq)) {
                seq.perm_reverse_complement();
           }
           (*poos) << seq.get_fasta();

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
