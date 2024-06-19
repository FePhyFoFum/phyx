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
#include "citations.h"


void print_help ();
std::string get_version_line ();

// write individual sequence to file named {sequence id}.fa
void write_sequence_to_file (Sequence& seq, std::vector<std::string>& outnames) {
    std::string fname = seq.get_id() + ".fa";
    // no spaces in file names!!!
    std::replace(fname.begin(), fname.end(), ' ', '_');
    outnames.push_back(fname);
    std::ofstream outf(fname);
    outf << ">" << seq.get_id() << std::endl;
    outf << seq.get_sequence() << std::endl;
    outf.close();
}

void print_help () {
    std::cout << std::endl;
    std::cout << "Split a multi-sequence alignment into separate files by taxon." << std::endl;
    std::cout << "This will take fasta, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << "Results are written in fasta format named '{sequence id}.fa'." << std::endl;
    std::cout << "Note: existing files will be overwritten." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxssplit [OPTIONS]... FILES" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input seq file, STDIN otherwise" << std::endl;
    std::cout << " -v, --verbose       make the output more verbose" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxssplit 1.3.1\n";
    vl += "Copyright (C) 2021 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"verbose", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool sfileset = false;
    bool verbose = false;
    char * seqf = nullptr;
    
    // keep track of the files produced
    std::vector<std::string> outnames;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:vhVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                sfileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
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
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    
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
    
    // this will be different for each taxon
//    if (outfileset) {
//        ofstr = new std::ofstream(outf);
//        poos = ofstr;
//    } else {
        poos = &std::cout;
//    }

    Sequence seq;
    std::string retstring;
    int num_taxa, num_char; // not used, but required by some reader functions

    int ft = test_seq_filetype_stream(*pios, retstring);

    // extra stuff to deal with possible interleaved nexus
    if (ft == 0) { // nexus
        bool interleave = false;
        get_nexus_dimensions(*pios, num_taxa, num_char, interleave);
        retstring = ""; // need to do this to let seqreader know we are mid-file
        if (!interleave) {
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                
                
                (*poos) << ">" << seq.get_id() << std::endl;
                (*poos) << seq.get_sequence() << std::endl;
            }
        } else {
            std::vector<Sequence> seqs = read_interleaved_nexus(*pios, num_taxa, num_char);
            for (const auto & sq : seqs) {
                seq = sq;
                
                
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
                
                
                (*poos) << ">" << seq.get_id() << std::endl;
                (*poos) << seq.get_sequence() << std::endl;
            }
        } else {
            // fasta, fastq, or simple phylip
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                write_sequence_to_file(seq, outnames);
            }
            // fasta has a trailing one
            if (ft == 2) {
                write_sequence_to_file(seq, outnames);
            }
        }
    }
    
    std::cerr << "Wrote " << outnames.size() << " files." << std::endl;
    if (verbose) {
        for (unsigned int i = 0; i < outnames.size(); i++) {
            std::cout << outnames[i] << std::endl;
        }
    }
    
    if (sfileset) {
        fstr->close();
        delete pios;
    }
    
    return EXIT_SUCCESS;
}
