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
#include "constants.h"

extern std::string PHYX_CITATION;


void print_help() {
    std::cout << "Taxon relabelling for alignments." << std::endl;
    std::cout << "This will take fasta, phylip or nexus file formats" << std::endl;
    std::cout << "Two ordered lists of taxa, -c (current) and -n (new) must be provided." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxrls [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -s, --seqf=FILE     input seq file, stdin otherwise" << std::endl;
    std::cout << " -c, --cnames=FILE   file containing current taxon labels (one per line)" << std::endl;
    std::cout << " -n, --nnames=FILE   file containing new taxon labels (one per line)" << std::endl;
    std::cout << " -o, --outf=FILE     output file, stout otherwise" << std::endl;
    std::cout << " -v, --verbose       make the output more verbose" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxrls 1.0\nCopyright (C) 2016-2019 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"cnames", required_argument, NULL, 'c'},
    {"nnames", required_argument, NULL, 'n'},
    {"outf", required_argument, NULL, 'o'},
    {"verbose", no_argument, NULL, 'v'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool sfileset = false;
    bool cfileset = false;
    bool nfileset = false;
    bool verbose = false;
    char * outf = NULL;
    char * seqf = NULL;
    std::string cnamef = "";
    std::string nnamef = "";
    
    while (1) {
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
                check_file_exists(cnamef.c_str());
                break;
            case 'n':
                nfileset = true;
                nnamef = strdup(optarg);
                check_file_exists(nnamef.c_str());
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
                print_error(argv[0], (char)c);
                exit(0);
        }
    }
    
    if (sfileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    std::istream * pios = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ofstream * ofstr = NULL;
    
    if (!nfileset | !cfileset) {
        std::cerr << "Error: must supply both name files (-c for current, -n for new). Exiting." << std::endl;
        exit(0);
    }
    
    if (sfileset == true) {
        fstr = new std::ifstream(seqf);
        pios = fstr;
    } else {
        pios = &std::cin;
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
            for (unsigned int i = 0; i < seqs.size(); i++) {
                seq = seqs[i];
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
            for (unsigned int i = 0; i < seqs.size(); i++) {
                seq = seqs[i];
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
    
    if (orig.size() > 0) {
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
