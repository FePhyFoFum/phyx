#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <cstring>
#include <getopt.h>

#include "aa2cdn.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"

void print_help() {
    std::cout << "Takes in an AA alignment and unaligned nucleotide file to generate a codon alignment." << std::endl;
    std::cout << "Taxa found in only 1 input file will be removed." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus inputs." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxaa2cdn [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -a, --aaseqf=FILE   input sequence file, stdin otherwise" << std::endl;
    std::cout << " -n, --nucseqf=FILE  input sequence file, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output fasta file, stout otherwise" << std::endl;
    std::cout << " -r, --rmlastcdn     removes last codon                " << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxaa2cdn 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"aaseqf", required_argument, NULL, 'a'},
    {"nucseqf", required_argument, NULL, 'n'},
    {"outf", required_argument, NULL, 'o'},
    {"rmlastcdn", no_argument, NULL, 'r'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    bool nucfileset = false;
    bool rm_last = false;
    char * aaseqf = NULL;
    char * nucseqf = NULL;
    char * outf = NULL;

    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "a:o:n:rhV", long_options, &oi);
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
    
    if (fileset && outfileset) {
        check_inout_streams_identical(aaseqf, outf);
    }
    if (nucfileset && outfileset) {
        check_inout_streams_identical(nucseqf, outf);
    }
    
    if (!fileset) {
        std::cerr << "You must specify an input amino acid sequence file. Exiting." << std::endl;
        exit(0);
    }
    if (!nucfileset) {
        std::cerr << "You must specify an input nucleotide sequence file. Exiting." << std::endl;
        exit(0);
    }
    
    std::ostream * poos = NULL;
    std::ofstream * ofstr = NULL;
    std::ifstream * fstr = NULL;
    std::istream * pios = NULL;
    std::ifstream * nucfstr = NULL;
    std::istream * nucpios = NULL;
    
    if (fileset == true) {
        fstr = new std::ifstream(aaseqf);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }
    if (nucfileset == true) {
        nucfstr = new std::ifstream(nucseqf);
        nucpios = nucfstr;
    } else {
        nucpios = &std::cin;
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
    
    Sequence aa_seq, nuc_seq;
    std::string retstring;
    
    
    // TODO: get rid of all of this map shit
    
    std::map<std::string, std::string> aa_sequences, nuc_sequences, codon_sequences;
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    while (read_next_seq_from_stream(*pios, ft, retstring, aa_seq)) {
        aa_sequences[aa_seq.get_id()] = aa_seq.get_sequence();
    }
    // fasta has a trailing one
    if (ft == 2) {
        aa_sequences[aa_seq.get_id()] = aa_seq.get_sequence();
    }
    
    // don't assume nuc and aa alignments are the same type
    ft = test_seq_filetype_stream(*nucpios, retstring);
    while (read_next_seq_from_stream(*nucpios, ft, retstring, nuc_seq)) {
        nuc_sequences[nuc_seq.get_id()] = nuc_seq.get_sequence();
    }
    // fasta has a trailing one
    if (ft == 2) {
        nuc_sequences[nuc_seq.get_id()] = nuc_seq.get_sequence();
    }
    
    AAtoCDN A2C;
    std::map<std::string, std::string>::iterator iter;
    codon_sequences = A2C.convert_to_codons(aa_sequences, nuc_sequences, rm_last);
    for (iter = codon_sequences.begin(); iter != codon_sequences.end(); iter++) {
        *poos << ">" << iter -> first << "\n" << iter -> second << std::endl;
    }
    
    
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
