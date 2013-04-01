/*
   Bare-bones sequence alignment resampling. Default is bootstrap, alternative is joackknife.
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <getopt.h>
#include <map>

using namespace std;

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "tree_reader.h"
#include "tree.h"
#include "cont_models.h"
#include "optimize_cont_models_nlopt.h"

void print_help(){
    cout << "Continuous character rate estimation with Brownian and OU" << endl;
    cout << "This will take fasta, phylip (and soon nexus) inputs." << endl;
    cout << endl;
    cout << "Usage: pxcontrates [OPTION]... " << endl;
    cout << endl;
    cout << " -c, --charf=FILE     input character file, stdin otherwise" << endl;
    cout << " -t, --treef=FILE     input tree file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxcontrates 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv2\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"char", required_argument, NULL, 'c'},
    {"tree", required_argument, NULL, 't'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    bool going = true;
    bool cfileset = false;
    bool tfileset = false;
    char * treef;
    char * charf;
    while (going) {
        int oi = -1;
        int c = getopt_long(argc, argv, "c:t:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'c':
                cfileset = true;
                charf = strdup(optarg);
                break;
            case 't':
                tfileset = true;
                treef = strdup(optarg);
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

    istream* pios;
    istream* poos;
    ifstream* cfstr;
    ifstream* tfstr;

    if (tfileset == true) {
        tfstr = new ifstream(treef);
        poos = tfstr;
    } 

    if (cfileset == true) {
        cfstr = new ifstream(charf);
        pios = cfstr;
    } 
    string retstring;
    int ft = test_char_filetype_stream(*pios, retstring);
    if(ft != 1 && ft != 2){
        cout << "only fasta and phylip (with spaces) supported so far" << endl;
        exit(0);
    }
    Sequence seq;
    vector<Sequence> seqs;
    map<string,int> seq_map;
    int y = 0;
    while (read_next_seq_char_from_stream(*pios,ft,retstring,seq)) {
        seqs.push_back(seq);
        seq_map[seq.get_id()] = y;
        seq.clear_cont_char();
        y++;
    }
    if(ft == 2){
        seqs.push_back(seq);
        seq_map[seq.get_id()] = y;
        seq.clear_cont_char();
    }
    //read trees
    TreeReader tr;
    vector<Tree *> trees;
    while(getline(*tfstr,retstring)){
        trees.push_back(tr.readTree(retstring));
    }

    //conduct analyses
    mat vcv;
    int t_ind = 0;
    int c_ind = 0;
    calc_vcv(trees[t_ind],vcv);
    int n = trees[t_ind]->getExternalNodeCount();
    rowvec x = rowvec(n);
    for (int i=0;i<n;i++){
        x(i) = seqs[seq_map[trees[t_ind]->getExternalNode(i)->getName()]].get_cont_char(c_ind);
    }
    vector<double> res = optimize_single_rate_bm_nlopt(x, vcv,true);
    cout << "state: " << res[0] <<  " rate: " << res[1] << " like: " << -res[2] << endl;;
    
    vector<double> res2 = optimize_single_rate_bm_ou_nlopt(x, vcv);
    cout << "state: " << res2[0] <<  " rate: " << res2[1] << " alpha: " << res2[2] <<  " like: " << -res2[3] << endl;;


    if (cfileset) {
        cfstr->close();
        delete pios;
    }
    if (tfileset) {
        tfstr->close();
        delete poos;
    }
}
