#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <getopt.h>

#include "utils.h"
#include "seq_reader.h"
#include "sequence.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Filter fastq files by mean quality." << std::endl;
    std::cout << "Data can be read from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxfqfilt [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -m, --mean=VALUE    mean value under which seqs are filtered" << std::endl;
    std::cout << " -s, --seqf=FILE     input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output sequence file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxfqfilt 1.3.1\n";
    vl += "Copyright (C) 2013-2024 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"mean", required_argument, nullptr, 'm'},
    {"seqf", required_argument, nullptr, 's'},
    {"outf", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    double meanfilt = 30;
    bool fileset = false;
    bool outfileset = false;
    char * seqf = nullptr;
    char * outf = nullptr;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "m:s:o:hVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'm':
                meanfilt = string_to_double(optarg, "-m");
                break;
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
    
    Sequence seq;
    std::string retstring;
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
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
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    if (ft != 3) {
        std::cerr << "Error: must be fastq input. Exiting." << std::endl;
        exit(1);
    }
    
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        double mean = seq.get_qualarr_mean();
        if (mean > meanfilt) {
            (*poos) << seq.get_fastq();
        }
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
