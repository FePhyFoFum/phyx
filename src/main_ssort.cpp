#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "seq_reader.h"
#include "sequence.h"
#include "seq_utils.h"
#include "utils.h"
#include "log.h"
#include "edlib.h"

void print_help() {
    cout << "Sort sequences by id or length" << endl;
    cout << "Can read from stdin or file, but output is fasta." << endl;
    cout << endl;
    cout << "Usage: pxssort [OPTION]... [FILE]..."<<endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise"<<endl;
    cout << " -b, --sortby        what to sort by: 1:id (default) 2:id rev" << endl;
    cout << "                                      3:length (<)   4:length (>)" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise"<<endl;
    cout << " -h, --help          display this help and exit"<<endl;
    cout << " -V, --version       display version and exit"<<endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" <<endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>"<<endl;
}

string versionline("pxssort 0.1\nCopyright (C) 2017 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"sortby", required_argument, NULL, 'b'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

struct SequenceIDListCompare{
    bool operator()(const Sequence & lhs, const Sequence & rhs){
      return lhs.get_id() < rhs.get_id();
    }
} SequenceIDListCompare;

struct SequenceRevIDListCompare{
    bool operator()(const Sequence & lhs, const Sequence & rhs){
      return lhs.get_id() > rhs.get_id();
    }
} SequenceRevIDListCompare;

struct SequenceLengthListCompare{
    bool operator()(const Sequence & lhs, const Sequence & rhs){
      return lhs.get_sequence().length() < rhs.get_sequence().length();
  }
} SequenceLengthListCompare;

struct SequenceRevLengthListCompare{
    bool operator()(const Sequence & lhs, const Sequence & rhs){
      return lhs.get_sequence().length() > rhs.get_sequence().length();
  }
} SequenceRevLengthListCompare;

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    int sortby = 1;
    
    char * seqf = NULL;
    char * outf = NULL;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:b:o:hgV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'b':
                sortby = atoi(strdup(optarg));
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
    vector<Sequence> seqs;
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios,retstring);
    while (read_next_seq_from_stream(*pios,ft,retstring,seq)) {
        seqs.push_back(seq);
    }
    if (ft == 2) {
        seqs.push_back(seq);
    }
    if (sortby == 1)
        sort(seqs.begin(),seqs.end(),SequenceIDListCompare);
    else if(sortby == 2)
        sort(seqs.begin(),seqs.end(),SequenceRevIDListCompare);
    else if(sortby == 3)
        sort(seqs.begin(),seqs.end(),SequenceLengthListCompare);
    else if(sortby == 4)
        sort(seqs.begin(),seqs.end(),SequenceRevLengthListCompare);
    for(unsigned int i=0;i<seqs.size();i++){
        (*poos) << seqs[i].get_fasta();
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
