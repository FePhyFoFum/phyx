#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <getopt.h>

#include "clean_seq.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"
#include "citations.h"


// TODO: throw out stop_codons: "TAG", "TAA", "TGA"
// TODO: read in partition file, edit and write out

void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Clean alignments by removing positions/taxa with too much ambiguous data." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << "Results are written in fasta format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxclsq [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -p, --prop=DOUBLE   proportion required to be present, default=0.5" << std::endl;
    std::cout << " -e, --empty         remove columns that are completely empty (- or ?)" << std::endl;
    std::cout << " -m, --min=INT       the minimum number of good characters required per site" << std::endl;
    std::cout << "                       - a min of 1 is equivalent to -e above" << std::endl;
    std::cout << " -t, --taxa          consider missing data per taxon (default: per site)" << std::endl;
    std::cout << " -c, --codon         examine sequences by codon rather than site" << std::endl;
    std::cout << "                       - requires all sequences be in frame and of correct length" << std::endl;
    std::cout << " -i, --info          report counts of missing data and exit" << std::endl;
    std::cout << "                       - combine with -t to get report by taxon (rather than site)" << std::endl;
    std::cout << "                       - combine with -c to use codons as units" << std::endl;
    std::cout << " -v, --verbose       more verbose output (i.e. if entire seqs are removed)" << std::endl;
    std::cout << " -o, --outf=FILE     output fasta file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxclsq 1.3.2\n";
    vl += "Copyright (C) 2015-2024 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"outf", required_argument, nullptr, 'o'},
    {"prop", required_argument, nullptr, 'p'},
    {"empty", no_argument, nullptr, 'e'},
    {"min", required_argument, nullptr, 'm'},
    {"taxa", required_argument, nullptr, 't'},
    {"codon", required_argument, nullptr, 'c'},
    {"info", required_argument, nullptr, 'i'},
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
    char * seqf = nullptr;
    char * outf = nullptr;
    double prop_required = 0.5;
    bool verbose = false;
    bool by_taxon = false;
    bool by_codon = false;
    bool count_only = false;
    bool remove_empty = false;
    int min_chars = 0;

    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:p:em:atcivhVC", long_options, &oi);
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
                prop_required = string_to_double(optarg, "-p");
                if (prop_required > 1.0 || prop_required < 0.0) {
                    std::cerr << "Error: proportion of required data present (-p) must be 0 <= p <= 1.0. Exiting."
                            << std::endl;
                    exit(0);
                }
                break;
            case 'e':
                remove_empty = true;
                break;
            case 'm':
                min_chars = string_to_int(optarg, "-m");;
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
    std::istream * pios = nullptr;
    std::ifstream * fstr = nullptr;
    
    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
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
    
    SequenceCleaner SC(pios, prop_required, remove_empty, min_chars, by_taxon,
            by_codon, count_only, verbose);
    
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
