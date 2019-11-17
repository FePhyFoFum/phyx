// this is essentially empty

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "log.h"

void print_help() {
    std::cout << "Basic Positive selection test" << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus inputs." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxkaks [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -i, --inputf=FILE      name of codon aligned fasta, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE        name of output file, stout otherwise" << std::endl;
    std::cout << "     --help             display this help and exit" << std::endl;
    std::cout << "     --version          display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxkaks 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"inputf", required_argument, NULL, 'i'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]){
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    char * outf = NULL;
    char * seqf = NULL;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "i:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'i':
                fileset = true;
                seqf = strdup(optarg);
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
                std::cout << "Error in input try -h" << std::endl;
                exit(0);
        }
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    std::istream * pios = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ofstream * ofstr = NULL;
    
    if (outfileset == true) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    if (fileset == true) {
        fstr = new ifstream(seqf);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }
    
    return EXIT_SUCCESS;
}
