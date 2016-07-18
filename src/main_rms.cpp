/*
 * main_rms.cpp
 *
 *  Created on: Jun 16, 2015
 *      Author: joe
 */


//g++ -std=c++11 rms.cpp main_rms.cpp utils.cpp superdouble.cpp sequence.cpp seq_reader.cpp seq_utils.cpp -o test
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <algorithm>

using namespace std;

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"

void print_help() {
    cout << "Removes unwanted sequences" << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxrms [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input nucleotide sequence file, stdin otherwise" << endl;
    cout << " -r, --rmf=FILE      input list of sequences to be removed each on a separate line" << endl;
    cout << " -o, --outf=FILE     output aa sequence file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxrms 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"rmf", required_argument, NULL, 'r'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
       
    bool fileset = false;
    bool rmfileset = false;
    bool outfileset = false;
    string seqf = "";
    string outf = "";
    string rmf = "";

    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:r:o:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'r':
                rmfileset = true;
                rmf = strdup(optarg);
                check_file_exists(rmf);
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
    ifstream* rstr;
    istream* rpios;
    
    if (fileset == true) {
        fstr = new ifstream(seqf);
        pios = fstr;
    } else {
        pios = &cin;
    }
    if (rmfileset == true) {
        rstr = new ifstream(rmf);
        rpios = rstr;
    } else {
        rpios = &cin;
    }
    
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    Sequence seq;
    string retstring;
    vector<string> to_remove;
    ifstream readline;
    string line;
    string seq_name;
    readline.open(rmf.c_str());
    if (readline.is_open()) {
        while (getline (readline, line)) {
            to_remove.push_back(line);
        }
    }

    int ft = test_seq_filetype_stream(*pios,retstring);
    //send sequences to be translated here
    while (read_next_seq_from_stream(*pios,ft,retstring,seq)) {
        seq_name = seq.get_id();
        if (find(to_remove.begin(), to_remove.end(), seq_name) != to_remove.end()) {
            
            // what is supposed to go here?
            
        } else {
            *poos << ">" << seq_name << "\n" << seq.get_sequence() << endl;
        }
    }
    //fasta has a trailing one
    if (ft == 2) {
        seq_name = seq.get_id();
        if (find(to_remove.begin(), to_remove.end(), seq_name) != to_remove.end()) {

        } else {
            *poos << ">" << seq_name << "\n" << seq.get_sequence() << endl;
        }
    }

    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
