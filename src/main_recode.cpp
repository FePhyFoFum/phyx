/*
 Bare-bones sequence recoding. RY-coding at first, but eventually codon-recoding.
 Codon-recoding will require genetic codes, and so knowledge of the taxon-specific codes.
 TODO: implement 'degen' coding.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "recode.h"
#include "log.h"

void print_help() {
    cout << "Nucleotide sequence recoding." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxrecode [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE      input sequence file, stdin otherwise" << endl;
    cout << " -r, --recode=STRING  string identifying recoding scheme (default: RY)" << endl;
    cout << "  Supported recodings (use any valid combination):" << endl;
    cout << "      R = A|G" << endl;
    cout << "      Y = C|T" << endl;
    cout << "      S = C|G" << endl;
    cout << "      W = A|T" << endl;
    cout << "      M = A|C" << endl;
    cout << "      K = G|T" << endl;
    cout << "      B = C|G|T" << endl;
    cout << "      D = A|G|T" << endl;
    cout << "      H = A|C|T" << endl;
    cout << "      V = A|C|G" << endl;
    cout << " -o, --outf=FILE      output sequence file, stout otherwise" << endl;
    cout << " -h, --help           display this help and exit" << endl;
    cout << " -V, --version        display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxrecode 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"recode", required_argument, NULL, 'r'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    string recodescheme = "";
    char * outf = NULL;
    char * seqf = NULL;
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
                recodescheme = strdup(optarg);
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
    
    // set default if arg not provided
    if (recodescheme == "") {
        recodescheme = "RY";
    }
    
    istream * pios = NULL;
    ostream * poos = NULL;
    ifstream * fstr = NULL;
    ofstream * ofstr = NULL;
    
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
    
    SequenceRecoder sr (recodescheme);
    
    Sequence seq;
    string retstring;
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        (*poos) << ">" << seq.get_id() << endl;
        (*poos) << sr.get_recoded_seq(seq.get_sequence()) << endl;
    }
// have to deal with last sequence outside while loop. fix this.
    if (ft == 2) {
        (*poos) << ">" << seq.get_id() << endl;
        (*poos) << sr.get_recoded_seq(seq.get_sequence()) << endl;
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
