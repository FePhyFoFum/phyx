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
#include "constants.h"

extern std::string PHYX_CITATION;


void print_help() {
    std::cout << "Generate a codon alignment from aligned amino acids and unaligned nucleotides." << std::endl;
    std::cout << "Taxa found in only 1 input file will be removed." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus inputs." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxaa2cdn [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -a, --aaseqf=FILE   input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -n, --nucseqf=FILE  input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output fasta file, STOUT otherwise" << std::endl;
    std::cout << " -r, --rmlastcdn     removes last codon (default: false)" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxaa2cdn 1.1\nCopyright (C) 2015-2020 FePhyFoFum\nLicense GPLv3\nWritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"aaseqf", required_argument, NULL, 'a'},
    {"nucseqf", required_argument, NULL, 'n'},
    {"outf", required_argument, NULL, 'o'},
    {"rmlastcdn", no_argument, NULL, 'r'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
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
        int c = getopt_long(argc, argv, "a:o:n:rhVC", long_options, &oi);
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
            case 'C':
                std::cout << PHYX_CITATION << std::endl;
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
        std::cerr << "Error: you must specify an input amino acid sequence file. Exiting." << std::endl;
        exit(0);
    }
    if (!nucfileset) {
        std::cerr << "Error: you must specify an input nucleotide sequence file. Exiting." << std::endl;
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
    // use general purpose reader
    std::vector<Sequence> nuc_seqs;
    std::vector<Sequence> aa_seqs;
    std::vector<Sequence> codon_seqs;
    std::string alphaName = "";
    
    // read in nucleotide seqs
    nuc_seqs = ingest_alignment(nucpios, alphaName);
    if (alphaName.compare("DNA") != 0) {
        std::cerr << "Error: incorrect alignment type provided. DNA was expected, but "
            << alphaName << " detected. Exiting." << std::endl;
        exit(0);
    }
    
    bool inFrame = true;
    unsigned int curN = 0;
    for (unsigned int i = 0; i < nuc_seqs.size(); i++) {
        curN = nuc_seqs[i].get_length();
        if (curN % 3 != 0) {
            std::cerr << "Error: nucleotide sequence length for '" << nuc_seqs[i].get_id()
                << "' is not a multiple of 3." << std::endl;
            inFrame = false;
        }
    }
    if (!inFrame) {
        std::cerr << "Error: nucleotide alignment does not appear to be in frame. Exiting" << std::endl;
        exit(0);
    }
    
    // and amino acid alignment
    aa_seqs = ingest_alignment(pios, alphaName);
    if (alphaName.compare("AA") != 0) {
        std::cerr << "Error: incorrect alignment type provided. Amino acids was expected, but "
            << alphaName << " detected. Exiting." << std::endl;
        exit(0);
    }
    
    AAtoCDN A2C(nuc_seqs, aa_seqs, rm_last);
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
