/*
 * main_upgma.cpp
 *
 *  Created on: Jun 10, 2015
 *      Author: joe
 */


// TODO: need to remove unnecessary includes
//g++ -std=c++11 upgma.cpp main_upgma.cpp utils.cpp superdouble.cpp sequence.cpp seq_reader.cpp seq_utils.cpp -o test
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <map>
#include <cstring>
#include <getopt.h>


using namespace std;

#include "upgma.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

void print_help() {
    cout << "Basic UPGMA Tree Maker." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxupgma [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output newick file, stout otherwise" << endl;
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
    
    ostream* poos;
    ofstream* ofstr;
    ifstream* fstr;
    istream* pios;
    if (fileset == true) {
        fstr = new ifstream(seqf);
        pios = fstr;
    } else {
        pios = &cin;
    }
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    UPGMA upgma(pios);
    *poos << upgma.get_newick() << endl;
    
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
