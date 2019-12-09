#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstring>
#include <getopt.h>
#include <sstream>

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "seq_info.h"
#include "log.h"
#include "constants.h"

extern std::string PHYX_CITATION;


void print_help() {
    std::cout << "Print sequence file summary." << std::endl;
    std::cout << "By default returns all properties. Alternatively choose 1 property." << std::endl;
    std::cout << "This will take fasta, phylip or nexus file formats" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxlssq [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -s, --seqf=FILE     input seq file, stdin otherwise" << std::endl;
    std::cout << " -i, --indiv         output stats for individual sequences" << std::endl;
    std::cout << " -n, --nseq          return the number of sequences" << std::endl;
    std::cout << " -c, --nchar         return the number of characters (only if aligned)" << std::endl;
    std::cout << "                        - for unaligned seqs, use with -i flag" << std::endl;
    std::cout << " -l, --labels        return all taxon labels (one per line)" << std::endl;
    std::cout << " -p, --prot          force interpret as protein (if inference fails)" << std::endl;
    std::cout << " -a, --aligned       return whether sequences are aligned (same length)" << std::endl;
    std::cout << " -f, --freqs         return character state frequencies" << std::endl;
    std::cout << " -m, --missing       return the proportion of missing (-,?) characters" << std::endl;
    std::cout << " -o, --outf=FILE     output stats file, stout otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxlssq 1.0\nCopyright (C) 2016-2019 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"indiv", no_argument, NULL, 'i'},
    {"nseq", no_argument, NULL, 'n'},
    {"nchar", no_argument, NULL, 'c'},
    {"labels", no_argument, NULL, 'l'},
    {"prot", no_argument, NULL, 'p'},
    {"aligned", no_argument, NULL, 'a'},
    {"freqs", no_argument, NULL, 'f'},
    {"missing", no_argument, NULL, 'm'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    bool indiv = false;
    bool force_protein = false; // i.e. if inference fails
    bool optionsset = false; // is true, do not return all properties
    bool get_labels = false;
    bool check_aligned = false;
    bool get_nseq = false;
    bool get_nchar = false;
    bool get_freqs = false;
    bool get_missing = false;
    char * outf = NULL;
    char * seqf = NULL;
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:inclpafmhVC", long_options, &oi);
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
            case 'i':
                indiv = true;
                break;
            case 'n':
                get_nseq = true;
                optionsset = true;
                break;
            case 'c':
                get_nchar = true;
                optionsset = true;
                break;
            case 'l':
                get_labels = true;
                optionsset = true;
                break;
            case 'p':
                force_protein = true;
                break;
            case 'a':
                check_aligned = true;
                optionsset = true;
                break;
            case 'f':
                get_freqs = true;
                optionsset = true;
                break;
            case 'm':
                get_missing = true;
                optionsset = true;
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
    std::ifstream * fstr = NULL;
    std::istream * pios = NULL;
    
    if ((get_labels + check_aligned + get_nseq + get_freqs + get_nchar + get_missing) > 1) {
        std::cerr << "Error: specify 1 property only (or leave blank to show all properties). Exiting." << std::endl;
        exit(0);
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
    if (outfileset == true) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    SeqInfo ls_Seq(pios, poos, indiv, force_protein);
    if (optionsset) {
        // get single property
        ls_Seq.get_property (get_labels, check_aligned, get_nseq,
            get_freqs, get_nchar, get_missing);
    } else {
        // the original behaviour
        ls_Seq.summarize();
    }
    
    if (fileset) {
        fstr->close();
        delete pios;
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    
    return EXIT_SUCCESS;
}

