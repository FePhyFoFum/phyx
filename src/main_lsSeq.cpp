//g++ -std=c++11 main_lsSeq.cpp ls_Seq.h ls_Seq.cpp utils.cpp utils.h superdouble.h superdouble.cpp sequence.h sequence.cpp seq_reader.h seq_reader.cpp seq_utils.h seq_utils.cpp -o test

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstring>
#include <getopt.h>
#include <sstream>


using namespace std;

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "ls_Seq.h"

void print_help() {
    cout << "Print sequence file summary" << endl;
    cout << "This will take fasta, phylip or nexus file formats" << endl;
    cout << endl;
    cout << "Usage: pxseqls [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input seq file, stdin otherwise" << endl;
    cout << " -a, --all=option    output stats of all sequences" << endl;
    cout << " -p, --prot=option   output stats for amino acids default is DNA" << endl;
    cout << " -o, --outf=FILE     output stats file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}
string versionline("pxlstr 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown and Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"all", no_argument, NULL, 'a'},
    {"prot", no_argument, NULL, 'p'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]){
	
	//cout << "Hello World\nMy name is Joe" << endl;
	bool outfileset = false;
    bool fileset = false;
    bool all = false;
    bool prot = false;
    char * outf;
    char * seqf;
    //string seqf = "";
    //string outf = "";
    
	while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:aphV", long_options, &oi);
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
            case 'a':
                all = true;
                break;
            case 'p':
                prot = true;
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
    
    Stats ls_Seq(pios, all, prot);
    
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

