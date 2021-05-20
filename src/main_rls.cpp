#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <cstring>
#include <getopt.h>

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "relabel.h"
#include "log.h"
#include "citations.h" // contains PHYX_CITATION


void print_help() {
    std::cout << "Taxon relabelling for alignments." << std::endl;
    std::cout << "This will take fasta, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << "Two ordered lists of taxa, -c (current) and -n (new) must be provided." << std::endl;
    std::cout << "Results are written in fasta format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxrls [OPTIONS]... FILES" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input seq file, STDIN otherwise" << std::endl;
    std::cout << " -c, --cnames=FILE   file containing current taxon labels (one per line)" << std::endl;
    std::cout << " -n, --nnames=FILE   file containing new taxon labels (one per line)" << std::endl;
    std::cout << " -v, --verbose       make the output more verbose" << std::endl;
    std::cout << " -o, --outf=FILE     output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxrls 1.2\nCopyright (C) 2016-2021 FePhyFoFum\nLicense GPLv3\nWritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"cnames", required_argument, nullptr, 'c'},
    {"nnames", required_argument, nullptr, 'n'},
    {"outf", required_argument, nullptr, 'o'},
    {"verbose", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool sfileset = false;
    bool cfileset = false;
    bool nfileset = false;
    bool verbose = false;
    char * outf = nullptr;
    char * seqf = nullptr;
    std::string cnamef;
    std::string nnamef;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:c:n:o:vhVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                sfileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'c':
                cfileset = true;
                cnamef = strdup(optarg);
                check_file_exists(cnamef);
                break;
            case 'n':
                nfileset = true;
                nnamef = strdup(optarg);
                check_file_exists(nnamef);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
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
                print_error(argv[0]);
                exit(0);
        }
    }
    
    if (sfileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
    if (!nfileset || !cfileset) {
        std::cerr << "Error: must supply both name files (-c for current, -n for new). Exiting." << std::endl;
        exit(0);
    }
    
    if (sfileset) {
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
    
    Relabel rl (cnamef, nnamef, verbose);
    
    std::set<std::string> orig = rl.get_names_to_replace();
    
    Sequence seq;
    std::string retstring;
    bool success = false;
    int num_taxa, num_char; // not used, but required by some reader functions
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    // extra stuff to deal with possible interleaved nexus
    if (ft == 0) {
        bool interleave = false;
        get_nexus_dimensions(*pios, num_taxa, num_char, interleave);
        retstring = ""; // need to do this to let seqreader know we are mid-file
        if (!interleave) {
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                std::string terp = seq.get_id();
                success = rl.relabel_sequence(seq);
                if (success) {
                    orig.erase(terp);
                }
                (*poos) << ">" << seq.get_id() << std::endl;
                (*poos) << seq.get_sequence() << std::endl;
            }
        } else {
            std::vector<Sequence> seqs = read_interleaved_nexus(*pios, num_taxa, num_char);
            for (const auto & sq : seqs) {
                seq = sq;
                std::string terp = seq.get_id();
                success = rl.relabel_sequence(seq);
                if (success) {
                    orig.erase(terp);
                }
                (*poos) << ">" << seq.get_id() << std::endl;
                (*poos) << seq.get_sequence() << std::endl;
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
                std::string terp = seq.get_id();
                success = rl.relabel_sequence(seq);
                if (success) {
                    orig.erase(terp);
                }
                (*poos) << ">" << seq.get_id() << std::endl;
                (*poos) << seq.get_sequence() << std::endl;
            }
        } else {
            // fasta, fastq, or simple phylip
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                std::string terp = seq.get_id();
                success = rl.relabel_sequence(seq);
                if (success) {
                    orig.erase(terp);
                }
                (*poos) << ">" << seq.get_id() << std::endl;
                (*poos) << seq.get_sequence() << std::endl;
            }
            // fasta has a trailing one
            if (ft == 2) {
                std::string terp = seq.get_id();
                success = rl.relabel_sequence(seq);
                if (success) {
                    orig.erase(terp);
                }
                (*poos) << ">" << seq.get_id() << std::endl;
                (*poos) << seq.get_sequence() << std::endl;
            }
        }
    }
    
    if (!orig.empty()) {
        if (verbose) {
            std::cerr << "The following names to match were not found in the alignment:" << std::endl;
            for (auto elem : orig) {
                std::cerr << elem << std::endl;
            }
        }
    }
    
    if (sfileset) {
        fstr->close();
        delete pios;
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
