#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cstring>
#include <getopt.h>

#ifdef OMP
#include <omp.h>
#endif

#include "utils.h"
#include "seq_reader.h"
#include "sequence.h"
#include "seq_utils.h"
#include "seq_models.h"
#include "pairwise_alignment.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Conduct Smith-Waterman analysis for all the seqs in a file." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << "Output is a list of the scores and distances (and the alignments if asked)." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxsw [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -a, --outalnf=FILE  output sequence file, won't output otherwise" << std::endl;
    std::cout << " -t, --seqtype=INT   sequence type, default=DNA (DNA=0,AA=1)" << std::endl;
    std::cout << " -m, --matrix=FILE   scoring matrix, default DNA=EDNAFULL, AA=BLOSUM62" << std::endl;
    std::cout << " -n, --nthreads=INT  number of threads (open mp), default=2" << std::endl;
    std::cout << " -v, --verbose       make the output more verbose, turns off parallel" << std::endl;
    std::cout << " -o, --outf=FILE     output score/distance file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxsw 1.3\n";
    vl += "Copyright (C) 2013-2021 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"outf", required_argument, nullptr, 'o'},
    {"outalnf", required_argument, nullptr, 'a'},
    {"seqtype", required_argument, nullptr, 't'},
    {"matrix", required_argument, nullptr, 'm'},
    {"nthreads", required_argument, nullptr, 'n'},
    {"verbose", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    bool outalnfileset = false;
    bool matrixfileset = false;
    char * seqf = nullptr;
    char * outf = nullptr;
    char * outaf = nullptr;
    char * matf = nullptr;
    int seqtype = 0; // DNA default, 1 = aa
    int num_threads = 2; // DNA default, 1 = aa
    bool verbose = false;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:a:t:m:n:vhVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'a':
                outalnfileset = true;
                outaf = strdup(optarg);
                break;
            case 't':
                seqtype = string_to_int(optarg, "-t");
                if (seqtype > 1) {
                    std::cerr << "Don't recognize seqtype " << seqtype
                            << ". Must be 0 (DNA) or 1 (AA)." << std::endl;
                    exit(0);
                }
                break;
            case 'm':
                matrixfileset = true;
                matf = strdup(optarg);
                check_file_exists(matf);
                break;
            case 'n':
                num_threads = string_to_int(optarg, "-n");
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
                std::cout << get_SW_citation() << std::endl;
                exit(0);
            default:
                print_error(*argv);
                exit(0);
        }
    }
    std::map<char, std::map<char, int> > sc_mat;
    if (matrixfileset) {
        read_scoring_matrix(matf, sc_mat);
    } else {
        if (seqtype == 0) {
            get_ednafull(sc_mat);
        } else { //aa
            
            // *** what is supposed to go here? ***
             
        }
    }
    Sequence seq;
    std::string retstring;
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    std::ofstream * afstr = nullptr;
    
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
    if (outalnfileset) {
        afstr = new std::ofstream(outaf);
    }

    std::string alphaName;
    std::vector<Sequence> seqs = ingest_alignment(pios, alphaName);

    // go all by all
    for (unsigned int i = 0; i < seqs.size(); i++) {
#ifdef OMP
        omp_set_num_threads(num_threads);
#endif
        #pragma omp parallel for
        for (unsigned int j = 0; j < seqs.size(); j++) {
            if (j > i) {
                std::string aln1;
                std::string aln2;
                double sc = sw(seqs[i], seqs[j], sc_mat, 0, aln1, aln2);
                #pragma omp critical
                {
                    (*poos) << seqs[i].get_id() << "\t" << seqs[j].get_id()
                            << "\t" << sc << std::endl;
                    if (verbose) {
                        std::cout << seqs[i].get_id() << "\t" << aln1 << "\n"
                                << seqs[j].get_id() << "\t" << aln2 << std::endl;
                    }
                    if (outalnfileset) {
                        (*afstr) << seqs[i].get_id() << "\t" << aln1 << "\n"
                                << seqs[j].get_id() << "\t" << aln2 << std::endl;
                    }
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
    if (outalnfileset) {
        afstr->close();
        delete afstr;
    }
    return EXIT_SUCCESS;
}
