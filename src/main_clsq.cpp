/*
 * main_clsq.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */



//compile alone: g++ -std=c++11 main_clsq.cpp clsq.cpp utils.cpp superdouble.cpp -o test
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

#include "clsq.h"
#include "utils.h"

void print_help() {
    cout << "Cleans up Sequences and Removes position with too much Ambiguous data." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxclsq [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output fasta file, stout otherwise" << endl;
    cout << " -p, --percent=DOUBLE  percent missing, default=50.0%" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxclsq 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
	{"percent", required_argument, NULL, 'p'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    bool fileset = false;
    bool outfileset = false;
    bool percentgiven = false;
    string seqf = "";
    string outf = ""; // not used at the moment
    double percent = 50.0;

    while (1) {
        int oi = -1;
        int curind = optind;
        int c = getopt_long(argc, argv,  "s:o:p:hV", long_options, &oi);
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
            case 'p':
            	percentgiven = true;
            	percent = atof(strdup(optarg));
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
    if (!fileset) {
        cout << "you must specify an input sequence file" << endl;
        exit(0);
    }

	clsq functions;
	ofstream ofstr(outf);
	string fasta;
	double MissingAllowed;
	map<string, string> sequences;
	map<string, string>::iterator iter;
	//HardCoded Test
	//fasta = ("TestFiles/clsq.test");
	fasta = seqf;
	MissingAllowed = percent;
	MissingAllowed = (MissingAllowed / 100.0);
	sequences = functions.FastaToOneLine(fasta, MissingAllowed);
    for(iter = sequences.begin(); iter != sequences.end(); iter++){
    	ofstr << ">" << iter -> first << "\n" << iter -> second << "\n";
    }

    return EXIT_SUCCESS;

}
