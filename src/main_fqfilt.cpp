#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <getopt.h>

#include "utils.h"
#include "seq_reader.h"
#include "sequence.h"
#include "log.h"

void print_help() {
    std::cout << "Filter fastq files by mean quality." << std::endl;
    std::cout << "Can read from stdin or file." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxfqfilt [OPTION]... [FILE]..."<<std::endl;
    std::cout << std::endl;
    std::cout << " -m, --mean=VALUE    mean value under which seqs are filtered" << std::endl;
    std::cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output sequence file, stout otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxfqfilt 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"mean", required_argument, NULL, 'm'},
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    double meanfilt = 30;
    bool fileset = false;
    bool outfileset = false;
    char * seqf = NULL;
    char * outf = NULL;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "m:s:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'm':
                meanfilt = string_to_float(optarg, "-m");
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
    
    Sequence seq;
    std::string retstring;
    
    std::istream * pios = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ofstream * ofstr = NULL;
    
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
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    if (ft != 3) {
        std::cerr << "Must be fastq input. Exiting." << std::endl;
        exit(1);
    }
    
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        double mean = seq.get_qualarr_mean();
        if (mean > meanfilt) {
            (*poos) << seq.get_fastq();
        }
    }
    
    if (fileset == true) {
        fstr->close();
        delete pios;
    }
    if (outfileset == true) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
