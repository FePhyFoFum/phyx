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

void print_help() {
    std::cout << "Translate DNA alignment to amino acids." << std::endl;
    std::cout << "NOTE: assumes sequences are in frame." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus inputs." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxtlate [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -s, --seqf=FILE     input nucleotide sequence file, stdin otherwise" << std::endl;
    std::cout << " -t, --table=STRING  which translation table to use. currently available:" << std::endl;
    std::cout << "                       'std' (standard, default)" << std::endl;
    std::cout << "                       'vmt' (vertebrate mtDNA)" << std::endl;
    std::cout << "                       'ivmt' (invertebrate mtDNA)" << std::endl;
    std::cout << "                       'ymt' (yeast mtDNA)" << std::endl;
    std::cout << " -o, --outf=FILE     output aa sequence file, stout otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxtlate 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");


static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"table", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    char * seqf = NULL;
    char * outf = NULL;
    std::string tab = "std";

    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:t:o:hV", long_options, &oi);
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
                std::cout << versionline << std::endl;
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
    std::ifstream * fstr = NULL;
    std::istream * pios = NULL;
    
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
    if (outfileset == true) {
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
    int ntax, nchar; // not used, but required by some readers
    std::string alphaName = ""; // will check first sequence type
    bool first = true;
    
    // extra stuff to deal with possible interleaved nexus
    if (ft == 0) {
        bool interleave = false;
        get_nexus_dimensions(*pios, ntax, nchar, interleave);
        retstring = ""; // need to do this to let seqreader know we are mid-file
        if (!interleave) {
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                if (first) {
                    alphaName = seq.get_alpha_name();
                    if (alphaName.compare("DNA") != 0) {
                        std::cerr << "Error: incorrect alignment type provided. DNA was expected, but "
                            << alphaName << " detected. Exiting." << std::endl;
                        exit(0);
                    }
                    first = false;
                }
                nuc_seq = seq.get_sequence();
                aa_seq = tl.translate(nuc_seq);
                (*poos) << ">" << seq.get_id() << "\n" << aa_seq << std::endl;
            }
        } else {
            std::vector<Sequence> seqs = read_interleaved_nexus(*pios, ntax, nchar);
            for (unsigned int i = 0; i < seqs.size(); i++) {
                seq = seqs[i];
                if (first) {
                    alphaName = seq.get_alpha_name();
                    if (alphaName.compare("DNA") != 0) {
                        std::cerr << "Error: incorrect alignment type provided. DNA was expected, but "
                            << alphaName << " detected. Exiting." << std::endl;
                        exit(0);
                    }
                    first = false;
                }
                nuc_seq = seq.get_sequence();
                aa_seq = tl.translate(nuc_seq);
                (*poos) << ">" << seq.get_id() << "\n" << aa_seq << std::endl;
            }
        }
    } else {
        bool complicated_phylip = false;
        // check if we are dealing with a complicated phylip format
        if (ft == 1) {
            get_phylip_dimensions(retstring, ntax, nchar);
            complicated_phylip = is_complicated_phylip(*pios, nchar);
        }
        if (complicated_phylip) {
            std::vector<Sequence> seqs = read_phylip(*pios, ntax, nchar);
            for (unsigned int i = 0; i < seqs.size(); i++) {
                seq = seqs[i];
                if (first) {
                    alphaName = seq.get_alpha_name();
                    if (alphaName.compare("DNA") != 0) {
                        std::cerr << "Error: incorrect alignment type provided. DNA was expected, but "
                            << alphaName << " detected. Exiting." << std::endl;
                        exit(0);
                    }
                    first = false;
                }
                nuc_seq = seq.get_sequence();
                aa_seq = tl.translate(nuc_seq);
                (*poos) << ">" << seq.get_id() << "\n" << aa_seq << std::endl;
            }
        } else {
            // fasta, fastq, or simple phylip
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                if (first) {
                    alphaName = seq.get_alpha_name();
                    if (alphaName.compare("DNA") != 0) {
                        std::cerr << "Error: incorrect alignment type provided. DNA was expected, but "
                            << alphaName << " detected. Exiting." << std::endl;
                        exit(0);
                    }
                    first = false;
                }
                nuc_seq = seq.get_sequence();
                aa_seq = tl.translate(nuc_seq);
                (*poos) << ">" << seq.get_id() << "\n" << aa_seq << std::endl;
            }
            // fasta has a trailing one
            if (ft == 2) {
                nuc_seq = seq.get_sequence();
                aa_seq = tl.translate(nuc_seq);
                (*poos) << ">" << seq.get_id() << "\n" << aa_seq << std::endl;
            }
        }
    }
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
