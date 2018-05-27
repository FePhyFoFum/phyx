
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "relabel.h"
#include "log.h"

void print_help() {
    cout << "Taxon relabelling for alignments." << endl;
    cout << "This will take fasta, phylip or nexus file formats" << endl;
    cout << "Two ordered lists of taxa, -c (current) and -n (new) must be provided." << endl;
    cout << endl;
    cout << "Usage: pxrls [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input seq file, stdin otherwise" << endl;
    cout << " -c, --cnames=FILE   file containing current taxon labels (one per line)" << endl;
    cout << " -n, --nnames=FILE   file containing new taxon labels (one per line)" << endl;
    cout << " -o, --outf=FILE     output file, stout otherwise" << endl;
    cout << " -v, --verbose       make the output more verbose" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxrls 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"cnames", required_argument, NULL, 'c'},
    {"nnames", required_argument, NULL, 'n'},
    {"outf", required_argument, NULL, 'o'},
    {"verbose", no_argument, NULL, 'v'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool sfileset = false;
    bool cfileset = false;
    bool nfileset = false;
    bool verbose = false;
    char * outf = NULL;
    char * seqf = NULL;
    string cnamef = "";
    string nnamef = "";
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:c:n:o:vhV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                sfileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'c':
                cfileset = true;
                cnamef = strdup(optarg);
                check_file_exists(cnamef.c_str());
                break;
            case 'n':
                nfileset = true;
                nnamef = strdup(optarg);
                check_file_exists(nnamef.c_str());
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
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
    
    if (sfileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    istream * pios = NULL;
    ostream * poos = NULL;
    ifstream * fstr = NULL;
    ofstream * ofstr = NULL;
    
    if (!nfileset | !cfileset) {
        cout << "Must supply both name files (-c for current, -n for new)." << endl;
        exit(0);
    }
    
    if (sfileset == true) {
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
    
    Relabel rl (cnamef, nnamef, verbose);
    
    set <string> orig = rl.get_names_to_replace();
    
    Sequence seq;
    string retstring;
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    bool success = false;
    
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        string terp = seq.get_id();
        success = rl.relabel_sequence(seq);
        if (success) {
            orig.erase(terp);
        }
        (*poos) << ">" << seq.get_id() << endl;
        (*poos) << seq.get_sequence() << endl;
    }
// have to deal with last sequence outside while loop. fix this.
    if (ft == 2) {
        string terp = seq.get_id();
        success = rl.relabel_sequence(seq);
        if (success) {
            orig.erase(terp);
        }
        (*poos) << ">" << seq.get_id() << endl;
        (*poos) << seq.get_sequence() << endl;
    }
    
    if (orig.size() > 0) {
        if (verbose) {
            cerr << "The following names to match were not found in the alignment:" << endl;
            for (auto elem : orig) {
                cerr << elem << endl;
            }
        }
    }
    
    if (sfileset) {
        fstr->close();
        delete pios;
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
