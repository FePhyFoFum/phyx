
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <getopt.h>

using namespace std;

#include "rate_model.h"
#include "tree.h"
#include "tree_reader.h"
#include "state_reconstructor_simple.h"
#include "seq_reader.h"
#include "sequence.h"
#include "seq_utils.h"
#include "utils.h"
#include "mcmc.h"

#include <armadillo>
using namespace arma;

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
    //read seqs
    vector<Sequence> seqs;
    vector<Sequence> sr_seqs;
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*spios,retstring);
    while(read_next_seq_from_stream(*spios,ft,retstring,seq)){
        (*poos) << seq.get_fasta();
	seqs.push_back(seq);
	Sequence tseq(seq.get_id(),"");
	sr_seqs.push_back(tseq);
    }
    if(ft == 2){
        (*poos) << seq.get_fasta();
	seqs.push_back(seq);
	Sequence tseq(seq.get_id(),"");
	sr_seqs.push_back(tseq);
    }
    

    //read trees 
    ft = test_tree_filetype_stream(*tpios, retstring);
    if(ft != 1){
	cerr << "this really only works with newick" << endl;
	exit(0);
    }going = true;
    Tree * tree;
    if(ft == 1){
	while(going){
	    tree = read_next_tree_from_stream_newick(*tpios,retstring,&going);
	    break;
	}
    }

    map<string,string> codon_dict;
    vector<string> codon_list;
    map<string,vector<int> > codon_pos;
    populate_codon_list(&codon_list);
    populate_map_codon_dict(&codon_dict);
    populate_map_codon_indices(&codon_pos);

    mat bf(61,61);
    mat K(61,61);
    mat w(61,61);
    mat inq(61,61);
    generate_bigpibf_K_w(&bf,&K,&w,codon_dict,codon_pos,codon_list);
    update_simple_goldman_yang_q(&inq,1.0,1.0,bf,K,w);
//    cout << inq << endl;

    RateModel rm(61);
    rm.set_Q(inq);
//    cout << rm.get_Q() << endl;
    int sites = (seqs[0].get_sequence().size()/3);
    StateReconstructorSimple sr(rm,sites);
    sr.set_tree(tree);
    cout << "there are " << sites << " sites" << endl;
    sm0_mcmc(375,10,tree,sr,rm,seqs,sr_seqs,codon_pos,bf,K,w,inq);
    
    sfstr->close();
    delete spios;
    tfstr->close();
    delete tpios;

    if(outfileset){
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}