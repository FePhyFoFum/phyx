/*
 * main_clsq.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "clsq.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"

void print_help() {
    cout << "Cleans alignments by removing positions with too much ambiguous data." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << "Results are written in fasta format." << endl;
    cout << endl;
    cout << "Usage: pxclsq [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE       input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE       output fasta file, stout otherwise" << endl;
    cout << " -p, --prop=DOUBLE     proportion required to be present, default=0.5" << endl;
    cout << " -a, --aminoacid       force interpret as protein (if inference fails)" << endl;
    cout << " -v, --verbose         more verbose output (i.e. if entire seqs are removed)" << endl;
    cout << " -h, --help            display this help and exit" << endl;
    cout << " -V, --version         display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxclsq 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"prop", required_argument, NULL, 'p'},
    {"aminoacid", required_argument, NULL, 'a'},
    {"verbose", no_argument, NULL, 'v'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    bool force_protein = false;
    char * seqf = NULL;
    char * outf = NULL;
    double proportion = 0.5;
    bool verbose = false;

    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:p:avhV", long_options, &oi);
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
                proportion = atof(strdup(optarg));
                break;
            case 'a':
                force_protein = true;
                break;
            case 'v':
                verbose = true;
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
    if (!fileset) {
        cout << "you must specify an input sequence file" << endl;
        exit(0);
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    ostream * poos = NULL;
    ofstream * ofstr = NULL;
    istream * pios = NULL;
    ifstream * fstr = NULL;
    
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
    
    SequenceCleaner toClean(pios, proportion, force_protein, verbose);
    
    // write sequences. currently only fasta format.
    toClean.write_seqs(poos);
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    if (fileset) {
        fstr->close();
        delete pios;
    }

    return EXIT_SUCCESS;
}
