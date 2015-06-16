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
#include <vector>
#include <sstream>
#include <iterator>
#include <cstring>
#include <getopt.h>


using namespace std;

#include "tlate.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

void print_help() {
    cout << "Basic Translate Program Uses Codon Table 1." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxtlate [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input nucleotide sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output aa sequence file, stout otherwise" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxupgma 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");


static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    bool fileset = false;
    bool outfileset = false;
    string seqf = "";
    string outf = "";

    while (1) {
        int oi = -1;
        int curind = optind;
        int c = getopt_long(argc, argv, "s:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
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
                cout << versionline << endl;
                exit(0);
            default:
                print_error(argv[0],(char)c);
                exit(0);
        }
    }

    // only taking files at the moment (not stdin)
    if (!fileset) {
        cout << "you must specify an input sequence file" << endl;
        exit(0);
    }
    //outfile prep
    ostream* poos;
    ofstream* ofstr;
    ifstream* fstr;
    istream* pios;
    if(fileset == true){
        fstr = new ifstream(seqf);
        pios = fstr;
    }else{
        pios = &cin;
    }
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    Sequence seq;
    TLATE functions;
    string retstring;
    string aa_seq;
    string nuc_seq;

    int ft = test_seq_filetype_stream(*pios,retstring);
    //send sequences to be translated here
    while(read_next_seq_from_stream(*pios,ft,retstring,seq)){
    	nuc_seq = seq.get_sequence();
    	aa_seq = functions.Translate(nuc_seq);
    	*poos << ">" << seq.get_id() << "\n" << aa_seq << endl;
    }
    //fasta has a trailing one
    if (ft == 2){
    	*poos << ">" << seq.get_id() << "\n" << aa_seq << endl;
    }

    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
