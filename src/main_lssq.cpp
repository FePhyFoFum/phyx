#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>
#include <sstream>

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "seq_info.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Print sequence file summary." << std::endl;
    std::cout << "By default returns all properties. Alternatively choose 1 property." << std::endl;
    std::cout << "This will take fasta, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxlssq [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input seq file, STDIN otherwise" << std::endl;
    std::cout << " -i, --indiv         output stats for individual sequences" << std::endl;
    std::cout << " -n, --nseq          return the number of sequences" << std::endl;
    std::cout << " -c, --nchar         return the number of characters (only if aligned)" << std::endl;
    std::cout << "                        - for unaligned seqs, use with -i flag" << std::endl;
    std::cout << " -l, --labels        return all taxon labels (one per line)" << std::endl;
    std::cout << " -p, --prot          force interpret as protein (if inference fails)" << std::endl;
    std::cout << " -a, --aligned       return whether sequences are aligned (same length)" << std::endl;
    std::cout << " -f, --freqs         return character state frequencies" << std::endl;
    std::cout << " -m, --missing       return the proportion of missing (-,?) characters" << std::endl;
    std::cout << " -o, --outf=FILE     output stats file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxlssq 1.3.2\n";
    vl += "Copyright (C) 2016-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph F. Walker, Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"outf", required_argument, nullptr, 'o'},
    {"indiv", no_argument, nullptr, 'i'},
    {"nseq", no_argument, nullptr, 'n'},
    {"nchar", no_argument, nullptr, 'c'},
    {"labels", no_argument, nullptr, 'l'},
    {"prot", no_argument, nullptr, 'p'},
    {"aligned", no_argument, nullptr, 'a'},
    {"freqs", no_argument, nullptr, 'f'},
    {"missing", no_argument, nullptr, 'm'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    bool indiv = false;
    bool force_protein = false; // i.e. if inference fails
    bool optionsset = false; // is true, do not return all properties
    int propcount = 0; //count how many properties are requested
    bool get_labels = false;
    bool check_aligned = false;
    bool get_nseq = false;
    bool get_nchar = false;
    bool get_freqs = false;
    bool get_missing = false;
    char * outf = nullptr;
    char * seqf = nullptr;
    
    while (true) {
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
                propcount++;
                break;
            case 'c':
                get_nchar = true;
                optionsset = true;
                propcount++;
                break;
            case 'l':
                get_labels = true;
                optionsset = true;
                propcount++;
                break;
            case 'p':
                force_protein = true;
                break;
            case 'a':
                check_aligned = true;
                optionsset = true;
                propcount++;
                break;
            case 'f':
                get_freqs = true;
                optionsset = true;
                propcount++;
                break;
            case 'm':
                get_missing = true;
                optionsset = true;
                propcount++;
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
    std::ifstream * fstr = nullptr;
    std::istream * pios = nullptr;
    
    if (propcount > 1) {
        std::cerr << "Error: specify 1 property only (or leave blank to show all properties). Exiting." << std::endl;
        exit(0);
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
    if (outfileset) {
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

