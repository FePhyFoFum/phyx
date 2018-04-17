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
    cout << " -n, --names=CSL     names sep by commas (NO SPACES!)" << endl;
    cout << " -f, --namesf=FILE   names in a file (each on a line)" << endl;
    cout << " -c, --comp          take the complement (i.e. move any taxa not in list)" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxrms 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"names",required_argument,NULL,'n'},
    {"namesf", required_argument, NULL, 'f'},
    {"comp", no_argument, NULL, 'c'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
       
    bool fileset = false;
    bool outfileset = false;
    bool namesset = false;
    bool namefileset = false;
    bool complement = false;
    char * namesc = NULL;
    char * namesfc = NULL;
    char * seqf = NULL;
    char * outf = NULL;
    string rmf = "";
    vector <string> names;

    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:n:f:co:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'n':
                namesset = true;
                namesc = strdup(optarg);
                break;
            case 'f':
                namefileset = true;
                namesfc = strdup(optarg);
                check_file_exists(namesfc);
                break;
            case 'c':
                complement = true;
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
    
    istream * pios = NULL;
    ostream * poos = NULL;
    ifstream * fstr = NULL;
    ofstream * ofstr = NULL;
    
    if (namesset == true) {
        vector <string> tokens2;
        string del2(",");
        tokens2.clear();
        tokenize(namesc, tokens2, del2);
        for (unsigned int j=0; j < tokens2.size(); j++) {
            trim_spaces(tokens2[j]);
            names.push_back(tokens2[j]);
        }
    } else if (namefileset == true) {
        ifstream nfstr(namesfc);
        string tline;
        while (getline(nfstr, tline)) {
            trim_spaces(tline);
            names.push_back(tline);
        }
        nfstr.close();
    } else {
        cerr << "you need to set the names of the taxa you want to remove (-n)" << endl;
        exit(0);
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
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    Sequence seq;
    string retstring;
    string seq_name;
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    if (!complement) {
        while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
            seq_name = seq.get_id();
            if (find(names.begin(), names.end(), seq_name) == names.end()) {
                *poos << ">" << seq_name << "\n" << seq.get_sequence() << endl;
            }
        }
        // fasta has a trailing one
        if (ft == 2) {
            seq_name = seq.get_id();
            if (find(names.begin(), names.end(), seq_name) == names.end()) {
                *poos << ">" << seq_name << "\n" << seq.get_sequence() << endl;
            }
        }
    } else {
        // keep the taxa passed in
        // complicating factor: not guaranteed to have any taxa left (i.e. empty output)
        while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
            seq_name = seq.get_id();
            if (find(names.begin(), names.end(), seq_name) != names.end()) {
                *poos << ">" << seq_name << "\n" << seq.get_sequence() << endl;
            }
        }
        // fasta has a trailing one
        if (ft == 2) {
            seq_name = seq.get_id();
            if (find(names.begin(), names.end(), seq_name) != names.end()) {
                *poos << ">" << seq_name << "\n" << seq.get_sequence() << endl;
            }
        }
    }

    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
