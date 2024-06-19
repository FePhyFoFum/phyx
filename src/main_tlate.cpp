#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <getopt.h>

#include "tlate.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Translate DNA alignment to amino acids." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << "NOTE: assumes sequences are in frame." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxtlate [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input nucleotide sequence file, STDIN otherwise" << std::endl;
    std::cout << " -t, --table=STRING  which translation table to use. currently available:" << std::endl;
    std::cout << "                       'std' (standard, default)" << std::endl;
    std::cout << "                       'vmt' (vertebrate mtDNA)" << std::endl;
    std::cout << "                       'ivmt' (invertebrate mtDNA)" << std::endl;
    std::cout << "                       'ymt' (yeast mtDNA)" << std::endl;
    std::cout << " -o, --outf=FILE     output aa sequence file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxtlate 1.3.1\n";
    vl += "Copyright (C) 2015-2024 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)";
    return vl;
}


static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"table", required_argument, nullptr, 't'},
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
    char * seqf = nullptr;
    char * outf = nullptr;
    std::string tab = "std";

    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:t:o:hVC", long_options, &oi);
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
                fileset = true;
                tab = strdup(optarg);
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
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    std::ostream * poos = nullptr;
    std::ofstream * ofstr = nullptr;
    std::ifstream * fstr = nullptr;
    std::istream * pios = nullptr;
    
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
    
    TLATE tl (tab);
    
    Sequence seq;
    std::string retstring;
    std::string aa_seq;
    std::string nuc_seq;
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    int num_taxa, num_char; // not used, but required by some readers
    std::string alphaName; // will check first sequence type
    bool first = true;
    
    // extra stuff to deal with possible interleaved nexus
    if (ft == 0) {
        bool interleave = false;
        get_nexus_dimensions(*pios, num_taxa, num_char, interleave);
        retstring = ""; // need to do this to let seqreader know we are mid-file
        if (!interleave) {
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                if (first) {
                    alphaName = seq.get_alpha_name();
                    if (alphaName != "DNA") {
                        std::cerr << "Error: incorrect alignment type provided. DNA was expected, but "
                            << alphaName << " detected. Exiting." << std::endl;
                        exit(0);
                    }
                    first = false;
                }
                nuc_seq = seq.get_sequence();
                aa_seq = tl.translate(nuc_seq);
                (*poos) << ">" << seq.get_id() << std::endl << aa_seq << std::endl;
            }
        } else {
            std::vector<Sequence> seqs = read_interleaved_nexus(*pios, num_taxa, num_char);
            for (const auto & sq : seqs) {
                seq = sq;
                if (first) {
                    alphaName = seq.get_alpha_name();
                    if (alphaName != "DNA") {
                        std::cerr << "Error: incorrect alignment type provided. DNA was expected, but "
                            << alphaName << " detected. Exiting." << std::endl;
                        exit(0);
                    }
                    first = false;
                }
                nuc_seq = seq.get_sequence();
                aa_seq = tl.translate(nuc_seq);
                (*poos) << ">" << seq.get_id() << std::endl << aa_seq << std::endl;
            }
        }
    } else {
        bool complicated_phylip = false;
        // check if we are dealing with a complicated phylip format
        if (ft == 1) {
            get_phylip_dimensions(retstring, num_taxa, num_char);
            complicated_phylip = is_complicated_phylip(*pios, num_char);
        }
        if (complicated_phylip) {
            std::vector<Sequence> seqs = read_phylip(*pios, num_taxa, num_char);
            for (const auto & sq : seqs) {
                seq = sq;
                if (first) {
                    alphaName = seq.get_alpha_name();
                    if (alphaName != "DNA") {
                        std::cerr << "Error: incorrect alignment type provided. DNA was expected, but "
                            << alphaName << " detected. Exiting." << std::endl;
                        exit(0);
                    }
                    first = false;
                }
                nuc_seq = seq.get_sequence();
                aa_seq = tl.translate(nuc_seq);
                (*poos) << ">" << seq.get_id() << std::endl << aa_seq << std::endl;
            }
        } else {
            // fasta, fastq, or simple phylip
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                if (first) {
                    alphaName = seq.get_alpha_name();
                    if (alphaName != "DNA") {
                        std::cerr << "Error: incorrect alignment type provided. DNA was expected, but "
                            << alphaName << " detected. Exiting." << std::endl;
                        exit(0);
                    }
                    first = false;
                }
                nuc_seq = seq.get_sequence();
                aa_seq = tl.translate(nuc_seq);
                (*poos) << ">" << seq.get_id() << std::endl << aa_seq << std::endl;
            }
            // fasta has a trailing one
            if (ft == 2) {
                nuc_seq = seq.get_sequence();
                aa_seq = tl.translate(nuc_seq);
                (*poos) << ">" << seq.get_id() << std::endl << aa_seq << std::endl;
            }
        }
    }
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
