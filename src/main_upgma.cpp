#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <getopt.h>

#include "upgma.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"
#include "citations.h" // contains PHYX_CITATION


void print_help() {
    std::cout << "Bare bones UPGMA tree generator" << std::endl;
    std::cout << "Currently only uses uncorrected p-distances." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxupgma [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output newick file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V,  --version      display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxupgma 1.2\nCopyright (C) 2015-2021 FePhyFoFum\nLicense GPLv3\nWritten by Joseph W. Brown, Joseph F. Walker");

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
    
    bool fileset = false;
    bool outfileset = false;
    char * seqf = nullptr;
    char * outf = nullptr;
    
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
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    std::ostream * poos = nullptr;
    std::ofstream * ofstr = nullptr;
    std::ifstream * fstr = nullptr;
    std::istream * pios = nullptr;
    
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
    
    UPGMA upgma(pios);
    (*poos) << upgma.get_newick() << std::endl;
    
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
