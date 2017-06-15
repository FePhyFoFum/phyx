
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
#include "comp_test.h"
#include "log.h"

void print_help() {
    cout << "Compositional homogeneity test" << endl;
    cout << "This will take fasta, phylip or nexus file formats" << endl;
    cout << endl;
    cout << "Usage: pxcomp [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input seq file, stdin otherwise" << endl;
    cout << " -p, --prot          force interpret as protein (if inference fails)" << endl;
    cout << " -o, --outf=FILE     output stats file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}
string versionline("pxcomp 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown and Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"prot", no_argument, NULL, 'p'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]){
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    char * outf = NULL;
    char * seqf = NULL;
    bool force_protein = false; // i.e. if inference fails
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:p:hV", long_options, &oi);
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
                force_protein = true;
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
    
    CompTest ct(pios, poos, force_protein);
    
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

