#include <iostream>
#include <string>
#include <fstream>
#include <vector>
//#include <sstream>
//#include <iterator>
#include <algorithm>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "sstat.h"
#include "log.h"

void print_help() {
    cout << "Calculates multinomial alignment test statistics." << endl;
    cout << "Currently only calculates the test statistic from Bollback (2002) MBE." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << endl;
    cout << "Usage: pxsstat [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input seq file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxsstat 0.1\nCopyright (C) 2017 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    char * outf = NULL;
    char * seqf = NULL;
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
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
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

    Sequence seq;
    vector<Sequence> seqs;
    string retstring;
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        seqs.push_back(seq);
    }
// have to deal with last sequence outside while loop. fix this.
    if (ft == 2) {
        seqs.push_back(seq);
    }
    
    //cout << "Read in " << seqs.size() << " sequences!" << endl;
    
    MultinomialSeqStat mm(seqs);
    (*poos) << mm.get_test_statistic() << endl;
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }

    return EXIT_SUCCESS;
}

