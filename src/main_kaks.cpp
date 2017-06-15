//g++ -std=c++11 main_kaks.cpp utils.cpp utils.h kaks.cpp kaks.h sequence.h sequence.cpp seq_utils.h seq_utils.cpp superdouble.cpp superdouble.h -o pxkaks

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "log.h"

void print_help() {
    cout << "Basic Positive selection test" << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxkaks [OPTION]... " << endl;
    cout << endl;
    cout << " -i, --inputf=FILE      name of codon aligned fasta, stdin otherwise" << endl;
    cout << " -o, --outf=FILE        name of output file, stout otherwise" << endl;
    cout << "     --help             display this help and exit" << endl;
    cout << "     --version          display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxkaks 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

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
                cout << versionline << endl;
                exit(0);
            default:
                cout << "Error in input try -h" << endl;
                exit(0);
        }
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    istream * pios = NULL;
    ostream * poos = NULL;
    ifstream * fstr = NULL;
    ofstream * ofstr = NULL;
    
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
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
    
    return EXIT_SUCCESS;
}
