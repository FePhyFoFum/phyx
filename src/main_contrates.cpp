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
    cout << "Continuous character rate estimation with Brownian and OU." << endl;
    cout << "This will take fasta, phylip (and soon nexus) inputs." << endl;
    cout << "Can read from stdin or file." << endl;
    cout << endl;
    cout << "Usage: pxcontrates [OPTION]... " << endl;
    cout << endl;
    cout << " -c, --charf=FILE     input character file, stdin otherwise" << endl;
    cout << " -t, --treef=FILE     input tree file, stdin otherwise" << endl;
    cout << " -a, --analysis=NUM   analysis type (0=anc[DEFAULT], 1=ratetest)" << endl;
    cout << " -o, --outf=FILE      output sequence file, stout otherwise" << endl;
    cout << "     --help           display this help and exit" << endl;
    cout << "     --version        display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxcontrates 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv2\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"char", required_argument, NULL, 'c'},
    {"tree", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
    {"analysis", required_argument, NULL, 'a'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    bool going = true;
    bool cfileset = false;
    bool tfileset = false;
    bool ofileset = false;
    
    char * treef;
    char * charf;
    char * outf;
    int analysis = 0;
    while (going) {
        int oi = -1;
        int c = getopt_long(argc, argv, "c:t:o:a:hV", long_options, &oi);
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
	        case 'o':
                ofileset = true;
                outf = strdup(optarg);
                break;
	        case 'a':
	        	if (optarg[0] == '1')
		            analysis = 1;
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
    }else{
        poos = &cin;
    }

    if (cfileset == true) {
        cfstr = new ifstream(charf);
        pios = cfstr;
    }else{
        cout << "you have to set a character file. Only a tree file can be read in through the stream;" << endl;
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
    int nchars = 0 ;
    while (read_next_seq_char_from_stream(*pios,ft,retstring,seq)) {
        seqs.push_back(seq);
        nchars = seq.get_num_cont_char();
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
    while(getline(*poos,retstring)){
        trees.push_back(tr.readTree(retstring));
    }

    //conduct analyses for each character
    for (int c =0;c < nchars;c++){
        cerr << "character: "<< c << endl;
	if (analysis == 0){
	    for(int i=0;i<trees[0]->getExternalNodeCount();i++){
		vector<Superdouble> tv (1);
		tv[0] = seqs[seq_map[trees[0]->getExternalNode(i)->getName()]].get_cont_char(c);
		trees[0]->getExternalNode(i)->assocDoubleVector("val",tv);
	    }
	    for(int i=0;i<trees[0]->getInternalNodeCount();i++){
		vector<Superdouble> tv (1);
		tv[0] = 0;
		trees[0]->getInternalNode(i)->assocDoubleVector("val",tv);
	    }
	    calc_square_change_anc_states(trees[0],0);
	    for(int i=0;i<trees[0]->getInternalNodeCount();i++){
		double tv = (*trees[0]->getInternalNode(i)->getDoubleVector("val"))[0];
		trees[0]->getInternalNode(i)->deleteDoubleVector("val");
		std::ostringstream s;
		s.precision(9);
		s << "[&value=" << tv << "]";
		trees[0]->getInternalNode(i)->setName(s.str());
	    }

	    for(int i=0;i<trees[0]->getExternalNodeCount();i++){
		double tv = (*trees[0]->getExternalNode(i)->getDoubleVector("val"))[0];
		trees[0]->getExternalNode(i)->deleteDoubleVector("val");
		std::ostringstream s;
		s.precision(9);
		s << fixed << trees[0]->getExternalNode(i)->getName() << "[&value=" << tv << "]";
		trees[0]->getExternalNode(i)->setName(s.str());
	    }
	    cout << "#nexus\nbegin trees;\ntree a =";
	    cout << trees[0]->getRoot()->getNewick(true);
	    cout << ";\nend;\n" << endl;
	}else if (analysis == 1){
	    mat vcv;
	    int t_ind = 0;//TODO: do this over trees
	    int c_ind = c;
	    calc_vcv(trees[t_ind],vcv);
	    int n = trees[t_ind]->getExternalNodeCount();
	    rowvec x = rowvec(n);
	    for (int i=0;i<n;i++){
		x(i) = seqs[seq_map[trees[t_ind]->getExternalNode(i)->getName()]].get_cont_char(c_ind);
	    }
	    vector<double> res = optimize_single_rate_bm_nlopt(x, vcv,true);
        double aic = (2*2)-(2*(-res[2]));
        double aicc = aic + ((2*2*(2+1))/(n-2-1));
	    cout << c << " BM "<< " state: " << res[0] <<  " rate: " << res[1] << " like: " << -res[2] << " aic: " << aic << " aicc: " << aicc <<  endl;
	    
	    vector<double> res2 = optimize_single_rate_bm_ou_nlopt(x, vcv);
        aic = (2*3)-(2*(-res2[3]));
        aicc = aic + ((2*3*(3+1))/(n-3-1));
	    cout << c << " OU "<< " state: " << res2[0] <<  " rate: " << res2[1] << " alpha: " << res2[2] <<  " like: " << -res2[3] << " aic: " << aic << " aicc: " << aicc << endl;
	}
    }

    if (cfileset) {
        cfstr->close();
        delete pios;
    }
    if (tfileset) {
        tfstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
