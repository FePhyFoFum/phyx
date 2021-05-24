#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "seq_reader.h"
#include "sequence.h"
#include "seq_utils.h"
#include "utils.h"
#include "log.h"
#include "edlib.h"
#include "citations.h" // contains PHYX_CITATION


void print_help (void);
std::string get_version_line (void);

void print_help () {
    std::cout << "Sort sequences by id or length." << std::endl;
    std::cout << "This will take fasta, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxssort [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -b, --sortby        what to sort by: 1:id (default) 2:id rev" << std::endl;
    std::cout << "                                      3:length (<)   4:length (>)" << std::endl;
    std::cout << " -o, --outf=FILE     output sequence file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxssort 1.2\n";
    vl += "Copyright (C) 2017-2021 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Stephen A. Smith (blackrim), Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"sortby", required_argument, nullptr, 'b'},
    {"outf", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

static struct SequenceIDListCompare {
    bool operator()(const Sequence& lhs, const Sequence& rhs) {
      return lhs.get_id() < rhs.get_id();
    }
} SequenceIDListCompare;

static struct SequenceRevIDListCompare {
    bool operator()(const Sequence& lhs, const Sequence& rhs) {
      return lhs.get_id() > rhs.get_id();
    }
} SequenceRevIDListCompare;

static struct SequenceLengthListCompare {
    bool operator()(const Sequence& lhs, const Sequence& rhs) {
      return lhs.get_sequence().length() < rhs.get_sequence().length();
  }
} SequenceLengthListCompare;

static struct SequenceRevLengthListCompare {
    bool operator()(const Sequence& lhs, const Sequence& rhs) {
      return lhs.get_sequence().length() > rhs.get_sequence().length();
  }
} SequenceRevLengthListCompare;

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    int sortby = 1;
    char * seqf = nullptr;
    char * outf = nullptr;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:b:o:hgVC", long_options,&oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'b':
                sortby = string_to_int(optarg, "-b");
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
            case 'V':
                std::cout << get_version_line() << std::endl;
                exit(0);
            case 'C':
                std::cout << PHYX_CITATION << std::endl;
                exit(0);
            default:
                print_error(argv[0]);
                exit(0);
        }
    }
    
    if (fileset&& outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
    if (fileset) {
        fstr = new std::ifstream(seqf);
        pios = fstr;
    } else {
        pios =&std::cin;
        if (!check_for_input_to_stream()) {
            print_help();
            exit(1);
        }
    }
    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos =&std::cout;
    }
    
    std::string alphaName;
    std::vector<Sequence> seqs = ingest_alignment(pios, alphaName);
    
    if (sortby == 1) {
        sort(seqs.begin(), seqs.end(), SequenceIDListCompare);
    } else if (sortby == 2) {
        sort(seqs.begin(), seqs.end(), SequenceRevIDListCompare);
    } else if (sortby == 3) {
        sort(seqs.begin(), seqs.end(), SequenceLengthListCompare);
    } else if (sortby == 4) {
        sort(seqs.begin(), seqs.end(), SequenceRevLengthListCompare);
    }
    for (auto & seq : seqs) {
        (*poos) << seq.get_fasta();
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
