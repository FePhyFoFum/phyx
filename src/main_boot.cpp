/*
 Bare-bones sequence alignment resampling. Default is bootstrap, alternative is jackknife.
 Conserved-partition bootstrap now implemented.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "seq_sample.h"
#include "log.h"
#include "constants.h" // contains PHYX_CITATION


void print_help() {
    std::cout << "Sequence alignment bootstrap or jackknife resampling." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxboot [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input sequence file, STDIN otherwise" << std::endl;
    std::cout << " -p, --partf=FILE    file listing empirical partitions: NAME = START-STOP[\\INTERVAL]" << std::endl;
    std::cout << " -f, --frac=DOUBLE   jackknife percentage, default bootstrap (i.e. -f 1.0)" << std::endl;
    std::cout << " -x, --seed=INT      random number seed, clock otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output sequence file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxboot 1.2\nCopyright (C) 2013-2021 FePhyFoFum\nLicense GPLv3\nWritten by Joseph W. Brown");

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"outf", required_argument, nullptr, 'o'},
    {"partf", required_argument, nullptr, 'p'},
    {"frac", required_argument, nullptr, 'f'},
    {"seed", required_argument, nullptr, 'x'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    bool partitioned = false;
    float jackfract = 0.0;
    char * outf = nullptr;
    char * seqf = nullptr;
    std::string partf = "";
    int seed = -1;
    while(true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:p:f:x:hVC", long_options, &oi);
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
            case 'p':
                partitioned = true;
                partf = strdup(optarg);
                break;
            case 'f':
                jackfract = string_to_float(optarg, "-f");
                if (jackfract < 0 || jackfract > 1) {
                    std::cerr << "Error: jackknife fraction must be 0 < x < 1. Exiting." << std::endl;
                    exit(0);
                }
                break;
            case 'x':
                seed = string_to_int(optarg, "-x");
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
    
    if (partitioned && (jackfract != 0.0)) {
        std::cerr << "Error: partitioned jackknife not implemented. Exiting." << std::endl;
        exit(0);
    }
    
    SequenceSampler ss(pios, seed, jackfract, partf);
    ss.write_resampled_seqs(poos);
    
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
