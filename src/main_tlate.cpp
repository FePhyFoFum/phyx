/*
 * main_tlate.cpp
 *
 *  Created on: Jun 16, 2015
 *      Author: joe
 */



//g++ -std=c++11 tlate.cpp main_tlate.cpp utils.cpp superdouble.cpp sequence.cpp seq_reader.cpp seq_utils.cpp -o test
#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "tlate.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"

void print_help() {
    cout << "Translate DNA alignment to amino acids." << endl;
    cout << "NOTE: assumes sequences are in frame." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxtlate [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input nucleotide sequence file, stdin otherwise" << endl;
    cout << " -t, --table=STRING  which translation table to use. currently available:" << endl;
    cout << "                       'std' (standard, default)" << endl;
    cout << "                       'vmt' (vertebrate mtDNA)" << endl;
    cout << "                       'ivmt' (invertebrate mtDNA)" << endl;
    cout << "                       'ymt' (yeast mtDNA)" << endl;
    cout << " -o, --outf=FILE     output aa sequence file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxtlate 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");


static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"table", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    char * seqf = NULL;
    char * outf = NULL;
    string tab = "std";

    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:t:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 't':
                fileset = true;
                tab = strdup(optarg);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
            case 'V':
                cout << versionline << endl;
                exit(0);
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    ostream * poos = NULL;
    ofstream * ofstr = NULL;
    ifstream * fstr = NULL;
    istream * pios = NULL;
    
    if (fileset == true) {
        fstr = new ifstream(seqf);
        pios = fstr;
    } else {
        pios = &cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    TLATE tl (tab);
    
    Sequence seq;
    string retstring;
    string aa_seq;
    string nuc_seq;

    int ft = test_seq_filetype_stream(*pios, retstring);
    // send sequences to be translated here
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        nuc_seq = seq.get_sequence();
        aa_seq = tl.translate(nuc_seq);
        (*poos) << ">" << seq.get_id() << "\n" << aa_seq << endl;
    }
    // fasta has a trailing one
    if (ft == 2) {
        nuc_seq = seq.get_sequence();
        aa_seq = tl.translate(nuc_seq);
        (*poos) << ">" << seq.get_id() << "\n" << aa_seq << endl;
    }

    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
