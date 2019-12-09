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
#include "constants.h"

extern std::string PHYX_CITATION;


// TODO: throw out stop_codons: "TAG", "TAA", "TGA"

void print_help() {
    std::cout << "Clean alignments by removing positions/taxa with too much ambiguous data." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus inputs." << std::endl;
    std::cout << "Results are written in fasta format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxclsq [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output fasta file, stout otherwise" << std::endl;
    std::cout << " -p, --prop=DOUBLE   proportion required to be present, default=0.5" << std::endl;
    std::cout << " -t, --taxa          consider missing data per taxon (default: per site)" << std::endl;
    std::cout << " -c, --codon         examine sequences by codon rather than site" << std::endl;
    std::cout << "                       - requires all sequences be in frame and of correct length" << std::endl;
    std::cout << " -i, --info          report counts of missing data and exit" << std::endl;
    std::cout << "                       - combine with -t to get report by taxon (rather than site)" << std::endl;
    std::cout << "                       - combine with -c to use codons as units" << std::endl;
    std::cout << " -v, --verbose       more verbose output (i.e. if entire seqs are removed)" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
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
    {"taxa", required_argument, NULL, 't'},
    {"codon", required_argument, NULL, 'c'},
    {"info", required_argument, NULL, 'i'},
    {"verbose", no_argument, NULL, 'v'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    char * seqf = NULL;
    char * outf = NULL;
    double prop_required = 0.5;
    bool verbose = false;
    bool by_taxon = false;
    bool by_codon = false;
    bool count_only = false;

    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:p:atcivhVC", long_options, &oi);
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
                prop_required = string_to_float(optarg, "-p");
                if (prop_required > 1.0 || prop_required < 0.0) {
                    std::cerr << "Error: proportion of required data present (-p) must be 0 <= p <= 1.0. Exiting." << std::endl;
                    exit(0);
                }
                break;
            case 't':
                by_taxon = true;
                break;
            case 'c':
                by_codon = true;
                break;
            case 'i':
                count_only = true;
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
            case 'C':
                std::cout << PHYX_CITATION << std::endl;
                exit(0);
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
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
    
    SequenceCleaner SC(pios, prop_required, by_taxon, by_codon, count_only, verbose);
    
    // write sequences. currently only fasta format.
    if (!count_only) {
        SC.write_seqs(poos);
    } else {
        SC.write_stats(poos);
    }
    
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
