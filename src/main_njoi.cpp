/*
 * main_NJOI.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: joe
 */




#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_procs() 1
    #define omp_get_max_threads() 1
#endif


// TODO: need to remove unnecessary includes
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <map>
#include <iterator>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "njoi.h"
#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"

//g++ -std=c++11 njoi.cpp main_njoi.cpp utils.cpp superdouble.cpp sequence.cpp seq_reader.cpp seq_utils.cpp -o test

void print_help() {
    cout << "Basic Neighbour-Joining Tree Maker." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxnjoi [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output newick file, stout otherwise" << endl;
    cout << " -n, --nthreads=INT  number of threads, default=1" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxnjoi 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");


static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"nthreads", required_argument, NULL, 'n'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

void printInfo () {
    int foo = omp_get_num_procs();
    int bar = omp_get_max_threads();
    cout << "numprocs = " << foo << "; threads = " << bar << endl;
}

int main(int argc, char * argv[]) {
    bool fileset = false;
    bool outfileset = false;
    int threads = 1;
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
            case 'n':
                threads = atoi(strdup(optarg));
                break;
            case 'h':
                print_help();
                printInfo(); // temp
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
//    if (!fileset) {
//        cout << "you must specify an input sequence file" << endl;
//        exit(0);
//    }
    
    //outfile prep
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
    
    NJOI nj(pios, threads);
    *poos << nj.get_newick() << endl;
    
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
