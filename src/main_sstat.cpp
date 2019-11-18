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

void print_help() {
    std::cout << "Calculates multinomial alignment test statistics." << std::endl;
    std::cout << "Currently only calculates the test statistic from Bollback (2002) MBE." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus inputs." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxsstat [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -s, --seqf=FILE     input seq file, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE     output sequence file, stout otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxsstat 0.1\nCopyright (C) 2017 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    char * outf = NULL;
    char * seqf = NULL;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:hV", long_options, &oi);
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
            default:
                print_error(argv[0], (char)c);
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
        fstr = new std::ifstream(seqf);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }

    Sequence seq;
    std::vector<Sequence> seqs;
    std::string retstring;
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        seqs.push_back(seq);
    }
// have to deal with last sequence outside while loop. fix this.
    if (ft == 2) {
        seqs.push_back(seq);
    }
    
    MultinomialSeqStat mm(seqs);
    (*poos) << mm.get_test_statistic() << std::endl;
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
