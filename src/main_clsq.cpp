#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <getopt.h>

#include "clsq.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"

void print_help() {
    std::cout << "Cleans alignments by removing positions with too much ambiguous data." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus inputs." << std::endl;
    std::cout << "Results are written in fasta format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxclsq [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -s, --seqf=FILE       input sequence file, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE       output fasta file, stout otherwise" << std::endl;
    std::cout << " -p, --prop=DOUBLE     proportion required to be present, default=0.5" << std::endl;
    std::cout << " -a, --aminoacid       force interpret as protein (if inference fails)" << std::endl;
    std::cout << " -v, --verbose         more verbose output (i.e. if entire seqs are removed)" << std::endl;
    std::cout << " -h, --help            display this help and exit" << std::endl;
    std::cout << " -V, --version         display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxclsq 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"prop", required_argument, NULL, 'p'},
    {"aminoacid", required_argument, NULL, 'a'},
    {"verbose", no_argument, NULL, 'v'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    bool force_protein = false;
    char * seqf = NULL;
    char * outf = NULL;
    double proportion = 0.5;
    bool verbose = false;

    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:p:avhV", long_options, &oi);
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
            case 'p':
                proportion = string_to_float(optarg, "-p");
                break;
            case 'a':
                force_protein = true;
                break;
            case 'v':
                verbose = true;
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
    if (!fileset) {
        std::cerr << "You must specify an input sequence file. Exiting." << std::endl;
        exit(0);
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    std::ostream * poos = NULL;
    std::ofstream * ofstr = NULL;
    std::istream * pios = NULL;
    std::ifstream * fstr = NULL;
    
    if (outfileset == true) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
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
    
    SequenceCleaner toClean(pios, proportion, force_protein, verbose);
    
    // write sequences. currently only fasta format.
    toClean.write_seqs(poos);
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    if (fileset) {
        fstr->close();
        delete pios;
    }

    return EXIT_SUCCESS;
}
