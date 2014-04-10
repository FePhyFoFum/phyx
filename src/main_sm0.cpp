
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <getopt.h>

using namespace std;

#include "seq_reader.h"
#include "sequence.h"
#include "seq_utils.h"
#include "utils.h"

void print_help(){
    cout << "Calculate Selection Model 0 with K and w with seqs and tree." << endl;
    cout << "Can read from stdin or file." << endl;
    cout << endl;
    cout << "Usage: pxsm0 [OPTION]... [FILE]..."<<endl;
    cout << endl; 
    cout << " -s, --seqf=FILE     input sequence file"<<endl;
    cout << " -t, --treef=FILE    input tree file" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise"<<endl;
    cout << "     --help          display this help and exit"<<endl;
    cout << "     --version       display version and exit"<<endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" <<endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>"<<endl;
}

string versionline("pxsm0 0.1\nCopyright (C) 2014 FePhyFoFum\nLicense GPLv2\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"treef", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]){
    bool going = true;
    bool sfileset = false;
    bool tfileset = false;
    bool outfileset = false;
    char * seqf;
    char * treef;
    char * outf;
    while(going){
        int oi = -1;
        int c = getopt_long(argc,argv,"s:t:o:hV",long_options,&oi);
        if (c == -1){
            break;
        }
        switch(c){
            case 's':
                sfileset = true;
                seqf = strdup(optarg);
                break;
	    case 't':
                tfileset = true;
                treef = strdup(optarg);
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
                print_error(argv[0],(char)c);
                exit(0);
        }
    }
    istream* spios;
    istream* tpios;
    ostream* poos;
    ifstream* sfstr;
    ifstream* tfstr;
    ofstream* ofstr; 
    if(sfileset == true){
        sfstr = new ifstream(seqf);
        spios = sfstr;
    }
    if(tfileset == true){
        tfstr = new ifstream(treef);
        tpios = tfstr;
    }
    if(tfileset == false || sfileset == false){
	cout << "you have to set the tree and seq files" << endl;
	exit(0);
    }
    if(outfileset == true){
        ofstr = new ofstream(outf);
        poos = ofstr;
    }else{
        poos = &cout;
    }
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*spios,retstring);
    while(read_next_seq_from_stream(*spios,ft,retstring,seq)){
        (*poos) << seq.get_fasta();
    }
    if(ft == 2){
        (*poos) << seq.get_fasta();
    }
    map<string,vector<int> > codon_pos;
    populate_map_codon_indices(&codon_pos);

    sfstr->close();
    delete spios;
    tfstr->close();
    delete tpios;

    if(outfileset){
        ofstr->close();
        delete poos;
    }
}
