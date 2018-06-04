
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>
#include <map>

using namespace std;

#include "string_node_object.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "tree_reader.h"
#include "tree.h"
#include "cont_models.h"
#include "optimize_cont_models_nlopt.h"
#include "log.h"

void print_help() {
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
    cout << " -h, --help           display this help and exit" << endl;
    cout << " -V, --version        display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxcontrates 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

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
    
    log_call(argc, argv);
    
    bool cfileset = false;
    bool tfileset = false;
    bool outfileset = false;
    
    char * treef = NULL;
    char * charf = NULL;
    char * outf = NULL;
    int analysis = 0;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "c:t:o:a:hV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'c':
                cfileset = true;
                charf = strdup(optarg);
                check_file_exists(charf);
                break;
            case 't':
                tfileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'a':
                if (optarg[0] == '1') {
                    analysis = 1;
                }
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

    istream * pios = NULL;
    istream * poos = NULL;
    ifstream * cfstr = NULL;
    ifstream * tfstr = NULL;

    ostream * poouts = NULL;
    ofstream * ofstr = NULL;
    

    if (tfileset == true) {
        tfstr = new ifstream(treef);
        poos = tfstr;
    } else {
        poos = &cin;
    }

    if (cfileset == true) {
        cfstr = new ifstream(charf);
        pios = cfstr;
    } else {
        cout << "you have to set a character file. Only a tree file can be read in through the stream;" << endl;
    }

    //out file
    //
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poouts = ofstr;
    } else {
        poouts = &cout;
    }
    //

    string retstring;
    int ft = test_char_filetype_stream(*pios, retstring);
    if (ft != 1 && ft != 2) {
        cout << "only fasta and phylip (with spaces) supported so far" << endl;
        exit(0);
    }
    Sequence seq;
    vector <Sequence> seqs;
    map <string, int> seq_map;
    int y = 0;
    int nchars = 0 ;
    while (read_next_seq_char_from_stream(*pios, ft, retstring, seq)) {
        seqs.push_back(seq);
        nchars = seq.get_num_cont_char();
        seq_map[seq.get_id()] = y;
        seq.clear_cont_char();
        y++;
    }
    
    if (ft == 2) {
        seqs.push_back(seq);
        seq_map[seq.get_id()] = y;
        seq.clear_cont_char();
    }
    //read trees
    TreeReader tr;
    vector<Tree *> trees;
    while (getline(*poos,retstring)) {
        trees.push_back(tr.readTree(retstring));
    }
    
    //conduct analyses for each character
    for (int c=0; c < nchars; c++) {
        cerr << "character: " << c << endl;
        if (analysis == 0) {
           // cout << "Input tree: " << getNewickString(trees[0]) << ";" << endl;
            if (c == 0) {
                (*poouts) << "#nexus" << endl << "begin trees;" << endl;
            }
            for (unsigned int x = 0; x < trees.size(); x++){
                for (int i=0; i < trees[x]->getExternalNodeCount(); i++) {
                    vector<Superdouble> tv (1);
                    tv[0] = seqs[seq_map[trees[x]->getExternalNode(i)->getName()]].get_cont_char(c);
                    trees[x]->getExternalNode(i)->assocDoubleVector("val",tv);
                }
                for (int i=0; i < trees[x]->getInternalNodeCount(); i++) {
                    vector<Superdouble> tv (1);
                    tv[0] = 0;
                    trees[x]->getInternalNode(i)->assocDoubleVector("val",tv);
                }
                calc_square_change_anc_states(trees[x],0); // second character dies here
                for (int i=0; i < trees[x]->getInternalNodeCount(); i++) {
                    double tv = (*trees[x]->getInternalNode(i)->getDoubleVector("val"))[0];
                    trees[x]->getInternalNode(i)->deleteDoubleVector("val");
                    std::ostringstream s;
                    s.precision(9);
                    s << fixed << tv;
                    StringNodeObject nob(s.str());
                    trees[x]->getInternalNode(i)->assocObject("value",nob);
                    //trees[x]->getInternalNode(i)->setName(s.str());
                }
                for (int i=0; i < trees[x]->getExternalNodeCount(); i++) {
                    double tv = (*trees[x]->getExternalNode(i)->getDoubleVector("val"))[0];
                    trees[x]->getExternalNode(i)->deleteDoubleVector("val");
                    std::ostringstream s;
                    s.precision(9);
                    s << fixed << tv;
                    StringNodeObject nob(s.str());
                    trees[x]->getExternalNode(i)->assocObject("value",nob);
                    //s << fixed << trees[x]->getExternalNode(i)->getName() << "[&value=" << tv << "]";
                    //trees[x]->getExternalNode(i)->setName(s.str());
                }
                (*poouts) << "tree tree" << c << " = ";
                (*poouts) << getNewickString(trees[x],"value") << endl;
            }
            if (c == (nchars - 1)) {
                (*poouts) << "end;\n" << endl;
            }
            // remove annotations
            remove_annotations(trees[0]);
            
        } else if (analysis == 1) {
            mat vcv;
            int t_ind = 0; // TODO: do this over trees
            int c_ind = c;
            calc_vcv(trees[t_ind],vcv);
            int n = trees[t_ind]->getExternalNodeCount();
            rowvec x = rowvec(n);
            for (int i=0; i < n; i++) {
                x(i) = seqs[seq_map[trees[t_ind]->getExternalNode(i)->getName()]].get_cont_char(c_ind);
            }
            vector<double> res = optimize_single_rate_bm_nlopt(x, vcv, true);
            double aic = (2*2)-(2*(-res[2]));
            double aicc = aic + ((2*2*(2+1))/(n-2-1));
            cout << c << " BM " << " state: " << res[0] <<  " rate: " << res[1]
                << " like: " << -res[2] << " aic: " << aic << " aicc: " << aicc <<  endl;

            vector<double> res2 = optimize_single_rate_bm_ou_nlopt(x, vcv);
            aic = (2*3)-(2*(-res2[3]));
            aicc = aic + ((2*3*(3+1))/(n-3-1));
            cout << c << " OU " << " state: " << res2[0] <<  " rate: "
                << res2[1] << " alpha: " << res2[2] <<  " like: " << -res2[3]
                << " aic: " << aic << " aicc: " << aicc << endl;
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
    if (outfileset) {
        ofstr->close();
        delete poouts;
    }
    return EXIT_SUCCESS;
}
