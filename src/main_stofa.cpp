/*
 * main_stofa.cpp
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "utils.h"
#include "seq_reader.h"
#include "sequence.h"

void print_help() {
    cout << "Convert seqfiles from nexus,phylip, fastq to fasta." << endl;
    cout << "Can read from stdin or file." << endl;
    cout << endl;
    cout << "Usage: pxstofa [OPTION]... [FILE]..." << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxstofa 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv2\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
	ofstream PhyxFile;
    PhyxFile.open ("Phyx.log", ios::out | ios::app );
	for (unsigned i = 0; i < argc;i++){
			PhyxFile << argv[i] << " ";
	}
	PhyxFile << endl;
	PhyxFile.close();
    bool fileset = false;
    bool outfileset = false;
    char * seqf;
    char * outf;
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
                cout << versionline << endl;
                exit(0);
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
    }
    Sequence seq;
    string retstring;
    istream* pios;
    ostream* poos;
    ifstream* fstr;
    ofstream* ofstr; 
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

    int ft = test_seq_filetype_stream(*pios,retstring);
    while (read_next_seq_from_stream(*pios,ft,retstring,seq)) {
        (*poos) << seq.get_fasta();
    }
    //fasta has a trailing one
    if (ft == 2) {
        (*poos) << seq.get_fasta();
    }
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

