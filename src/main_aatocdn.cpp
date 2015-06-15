/*
 * main_aatocdn.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */


//To compile alone: g++ -std=c++11 main_aatocdn.cpp aatocdn.cpp utils.cpp superdouble.cpp -o test
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

#include "aatocdn.h"
#include "utils.h"

void print_help() {
    cout << "Takes in AA alignment and unaligned Nucleotide file to give Codon Alignment." << endl;
    cout << "This will get rid of any sequences found in either only the Nucleotide or the Amino Acid Alignment" << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxaatocdn [OPTION]... " << endl;
    cout << endl;
    cout << " -a, --aaseqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -n, --nucseqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output fasta file, stout otherwise" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxclsq 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"aaseqf", required_argument, NULL, 'a'},
	{"nucseqf", required_argument, NULL, 'n'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    bool fileset = false;
    bool outfileset = false;
    bool nucfile = false;
    string aaseqf = "";
    string nucseqf = "";
    string outf = "";

    while (1) {
        int oi = -1;
        int curind = optind;
        int c = getopt_long(argc, argv,  "a:o:n:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'a':
                fileset = true;
                aaseqf = strdup(optarg);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'n':
            	nucfile = true;
                nucseqf = strdup(optarg);
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
    ostream* poos;
    ofstream* ofstr;
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }

	aa_to_cdn functions;
	string AAFasta, NucFasta;
	map<string, string> AASequences, NucSequences, CodonSequences;
	map<string, string>::iterator iter;
	AAFasta = ("TestFiles/AA.fa");
	NucFasta = ("TestFiles/Codon.fa");
	AASequences = functions.FastaToOneLine(AAFasta);
	NucSequences = functions.FastaToOneLine(NucFasta);
	CodonSequences = functions.ChangeToCodon(AASequences, NucSequences);
    for (iter = CodonSequences.begin(); iter != CodonSequences.end(); iter++){
    	*poos << ">" << iter -> first << "\n" << iter -> second << endl;
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;

}


