#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <getopt.h>

#include "aa2cdn.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Generate a codon alignment from aligned amino acids and unaligned nucleotides." << std::endl;
    std::cout << "Taxa found in only 1 input file will be removed." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus inputs." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxaa2cdn [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -a, --aaseqf=FILE   input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -n, --nucseqf=FILE  input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -r, --rmlastcdn     remove last codon from *all* nuc sequences (default: false)" << std::endl;
    std::cout << " -s, --stopremove    remove stop codon from nuc sequences if present (default: false)" << std::endl;
    std::cout << " -o, --outf=FILE     output fasta file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxaa2cdn 1.3.2\n";
    vl += "Copyright (C) 2015-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"aaseqf", required_argument, nullptr, 'a'},
    {"nucseqf", required_argument, nullptr, 'n'},
    {"outf", required_argument, nullptr, 'o'},
    {"rmlastcdn", no_argument, nullptr, 'r'},
    {"stopremove", no_argument, nullptr, 's'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    bool nucfileset = false;
    bool rm_last = false;
    bool rm_stop = false;
    char * aaseqf = nullptr;
    char * nucseqf = nullptr;
    char * outf = nullptr;


    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "a:o:n:rshVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'a':
                fileset = true;
                aaseqf = strdup(optarg);
                check_file_exists(aaseqf);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'n':
                nucfileset = true;
                nucseqf = strdup(optarg);
                check_file_exists(nucseqf);
                break;
            case 'r':
                rm_last = true;
                break;
            case 's':
                rm_stop = true;
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
        check_inout_streams_identical(aaseqf, outf);
    }
    if (nucfileset && outfileset) {
        check_inout_streams_identical(nucseqf, outf);
    }
    
    std::ostream * poos = nullptr;
    std::ofstream * ofstr = nullptr;
    std::ifstream * fstr = nullptr;
    std::istream * pios = nullptr;
    std::ifstream * nucfstr = nullptr;
    std::istream * nucpios = nullptr;
    
    if (fileset) {
        fstr = new std::ifstream(aaseqf);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (!check_for_input_to_stream()) {
            print_help();
            exit(1);
        }
    }
    if (nucfileset) {
        nucfstr = new std::ifstream(nucseqf);
        nucpios = nucfstr;
    } else {
        nucpios = &std::cin;
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
    
    if (!fileset) {
        std::cerr << "Error: you must specify an input amino acid sequence file. Exiting." << std::endl;
        exit(0);
    }
    if (!nucfileset) {
        std::cerr << "Error: you must specify an input nucleotide sequence file. Exiting." << std::endl;
        exit(0);
    }
    
    // use general purpose reader
    std::vector<Sequence> nuc_seqs;
    std::vector<Sequence> aa_seqs;
    std::string alphaName;
    
    // read in nucleotide seqs
    nuc_seqs = ingest_alignment(nucpios, alphaName);
    if (alphaName != "DNA") {
        std::cerr << "Error: incorrect alignment type provided. DNA was expected, but "
            << alphaName << " detected. Exiting." << std::endl;
        exit(0);
    }
    
    bool inFrame = true;
    for (auto & nuc_seq : nuc_seqs) {
        unsigned int  curN = nuc_seq.get_length();
        if (curN % 3 != 0) {
            std::cerr << "Error: nucleotide sequence length for '" << nuc_seq.get_id()
                << "' is not a multiple of 3." << std::endl;
            inFrame = false;
        }
    }
    if (!inFrame) {
        std::cerr << "Error: nucleotide alignment does not appear to be in frame. Exiting." << std::endl;
        exit(0);
    }
    
    // and amino acid alignment
    aa_seqs = ingest_alignment(pios, alphaName);
    if (alphaName != "AA") {
        std::cerr << "Error: incorrect alignment type provided. Amino acids was expected, but "
            << alphaName << " detected. Exiting." << std::endl;
        exit(0);
    }
    
    if (rm_last && rm_stop) {
        std::cerr << "Error: you may set -r or -s, but not both. Exiting." << std::endl;
        exit(0);
    }
    
    AAtoCDN A2C(nuc_seqs, aa_seqs, rm_last, rm_stop);
    A2C.write_codon_alignment(poos);
    
    if (fileset) {
        fstr->close();
        delete pios;
    }
    if (nucfileset) {
        nucfstr->close();
        delete nucpios;
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
