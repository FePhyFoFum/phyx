#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <getopt.h>

#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "sstat.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Calculates multinomial alignment test statistics." << std::endl;
    std::cout << "Currently only calculates the test statistic from Bollback (2002) MBE." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxsstat [OPTIONS]... " << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input seq file, STDIN otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output sequence file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxsstat 1.3.2\n";
    vl += "Copyright (C) 2017-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"outf", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    char * outf = nullptr;
    char * seqf = nullptr;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:hVC", long_options, &oi);
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
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
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

    std::string alphaName;
    std::vector<Sequence> seqs = ingest_alignment(pios, alphaName);
    
    MultinomialSeqStat mm(seqs);
    (*poos) << mm.get_test_statistic() << std::endl;
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
