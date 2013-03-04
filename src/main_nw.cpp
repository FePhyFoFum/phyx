/*
 * main_mrca.cpp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <getopt.h>

using namespace std;

#include "utils.h"
#include "seq_reader.h"
#include "sequence.h"
#include "seq_utils.h"

void print_help(){
    cout << "Conduct Needleman-Wunsch analysis for all the seqs in a file." << endl;
	cout << "This will take fasta, fastq, phylip, and nexus inputs. And " << endl;
	cout << "will output a list of the scores and distances (and the alignments"<< endl;
	cout << "if asked)." << endl;
    cout << "Can read from stdin or file." << endl;
    cout << endl;
    cout << "Usage: pxnw [OPTION]... [FILE]..."<<endl;
    cout << endl; 
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise"<<endl;
    cout << " -o, --outf=FILE     output score/distance file, stout otherwise"<<endl;
    cout << " -a, --outalnf=FILE  output sequence file, won't output otherwise"<<endl;
    cout << " -t, --seqtype=INT   sequence type, default=DNA (DNA=0,AA=1)"<<endl;
    cout << " -m, --matrix=FILE   scoring matrix, default DNA=EDNAFULL, AA=BLOSUM62"<<endl;
    cout << "     --help          display this help and exit"<<endl;
    cout << "     --version       display version and exit"<<endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" <<endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>"<<endl;
}

string versionline("pxnw 0.1\nCopyright (C) 2013 FePhyFoFum\nLiscence GPLv2\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]){
    bool going = true;
    bool fileset = false;
    bool outfileset = false;
	bool outalnfileset = false;
	bool matrixfileset = false;
    char * seqf;
    char * outf;
	char * outaf;
	char * matf;
	int seqtype = 0;//DNA default, 1 = aa
    while(going){
        int oi = -1;
        int c = getopt_long(argc,argv,"s:o:a:t:m:hV",long_options,&oi);
        if (c == -1){
            break;
        }
        switch(c){
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'a':
                outalnfileset = true;
                outaf = strdup(optarg);
                break;
            case 't':
                seqtype = atoi(strdup(optarg));
				if(seqtype > 1){
					cout << "Don't recognize seqtype " << seqtype << ". Must be 0 (DNA) or 1 (AA)." << endl;
					exit(0);
				}
                break;
            case 'm':
                matrixfileset = true;
                matf = strdup(optarg);
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
    vector<Sequence> seqs;
    Sequence seq;
    string retstring;
    istream* pios;
    ostream* poos;
    ifstream* fstr;
    ofstream* ofstr; 
    if(fileset == true){
        fstr = new ifstream(seqf);
        pios = fstr;
    }else{
        pios = &cin;
    }
    if(outfileset == true){
        ofstr = new ofstream(outf);
        poos = ofstr;
    }else{
        poos = &cout;
    }

	int ft = test_seq_filetype_stream(*pios,retstring);
	while(read_next_seq_from_stream(*pios,ft,retstring,seq)){
	    seqs.push_back(seq);
	}
	//fasta has a trailing one
	if (ft == 2){
	    seqs.push_back(seq);
	}
    write_phylip_alignment(seqs,poos);
    if(fileset){
        fstr->close();
        delete pios;
    }if(outfileset){
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}

