/*
 * main_NJ.cpp
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


#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "nj.h"
#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "log.h"

//g++ -std=c++11 nj.cpp main_nj.cpp utils.cpp superdouble.cpp sequence.cpp seq_reader.cpp seq_utils.cpp log.cpp -o test

void print_help() {
    cout << "Basic Neighbour-Joining Tree Maker." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxnj [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output newick file, stout otherwise" << endl;
    cout << " -n, --nthreads=INT  number of threads, default=1" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxnj 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");


static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"nthreads", required_argument, NULL, 'n'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

// temp function to play around with multithreading
void printInfo () {
    int foo = omp_get_num_procs();
    int bar = omp_get_max_threads();
    cout << "numprocs = " << foo << "; threads = " << bar << endl;
}

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    int threads = 1;
    char * seqf = NULL;
    char * outf = NULL;

    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:n:hV", long_options, &oi);
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
            case 'n':
                threads = atoi(strdup(optarg));
                break;
            case 'h':
                print_help();
                //printInfo(); // temp
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
