#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <time.h>
#include <stdlib.h>
#include <getopt.h>

using namespace std;

#include "seq_reader.h"
#include "sequence.h"
#include "seq_utils.h"
#include "utils.h"
#include "log.h"
#include "edlib.h"

void print_help() {
    cout << "Reverse complement sequences from nexus, phylip, or fastq to fasta." << endl;
    cout << "Can read from stdin or file, but output is fasta." << endl;
    cout << endl;
    cout << "Usage: pxrevcomp [OPTION]... [FILE]..."<<endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise"<<endl;
    cout << " -i, --ids=IDS       a comma sep list of ids to flip (NO SPACES!)" << endl;
    cout << " -g, --guess         EXPERIMENTAL: guess whether there are seqs that need to be " << endl;
    cout << "                       rev comp. uses edlib library on first seq" << endl;
    cout << " -p, --pguess        EXPERIMENTAL: progressively guess " << endl;
    cout << " -m, --sguess        EXPERIMENTAL: sampled guess " << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise"<<endl;
    cout << " -h, --help          display this help and exit"<<endl;
    cout << " -V, --version       display version and exit"<<endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" <<endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>"<<endl;
}

string versionline("pxrevcomp 0.11\nCopyright (C) 2017 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"ids", required_argument, NULL, 'i'},
    {"guess", no_argument, NULL, 'g'},
    {"pguess", no_argument, NULL, 'p'},
    {"sguess", no_argument, NULL, 'm'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

bool reverse_it_or_not(vector<Sequence> & seqs, Sequence comp_seq){
    int best_distance = 10000000;
    int best_dis_rev = 100000000;
    string comp = comp_seq.get_sequence();
    string revcomp = comp_seq.reverse_complement();
    for (unsigned int i=0;i<seqs.size();i++){
        EdlibAlignResult result = edlibAlign(comp.c_str(), comp.length(), seqs[i].get_sequence().c_str(), 
                seqs[i].get_sequence().length(), edlibNewAlignConfig(best_distance, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
        if (result.editDistance < 0)
            continue;
        if (result.editDistance < best_distance)
            best_distance = result.editDistance;
        edlibFreeAlignResult(result);
    }
    for (unsigned int i=0;i<seqs.size();i++){
        EdlibAlignResult result = edlibAlign(revcomp.c_str(), revcomp.length(), seqs[i].get_sequence().c_str(), 
                seqs[i].get_sequence().length(), edlibNewAlignConfig(best_distance, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
        if (result.editDistance < 0)
            continue;
        if (result.editDistance < best_dis_rev)
            best_dis_rev = result.editDistance;
        edlibFreeAlignResult(result);
    }
    if (best_dis_rev < best_distance)
        return true;
    return false;

}

int main(int argc, char * argv[]) {
    srand (time(NULL));
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    bool idsset = false;
    vector<string> ids;
    
    bool guess = false;
    bool pguess = false;
    bool sguess = false;
    double sguess_samplenum = 0.2; // 10% of them will be used for revcomp
    char * seqf = NULL;
    char * outf = NULL;
    char * idssc = NULL;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:i:o:mgphV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'i':
                idsset = true;
                idssc = strdup(optarg);
                break;
            case 'g':
                guess = true;
                break;
            case 'p':
                guess = true;
                pguess = true;
                break;
            case 'm':
                guess = true;
                sguess = true;
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
    
    if (idsset == true){
        vector<string> tokens2;
        tokenize(idssc, tokens2, ",");
        for (unsigned int j=0; j < tokens2.size(); j++) {
            trim_spaces(tokens2[j]);
            ids.push_back(tokens2[j]);
        }
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
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios,retstring);
    if (guess == false){
        while (read_next_seq_from_stream(*pios,ft,retstring,seq)) {
            if(idsset == false || count(ids.begin(),ids.end(),seq.get_id())==1){
                seq.perm_reverse_complement();
            }
            (*poos) << seq.get_fasta();
        }
        if (ft == 2) {
            if(idsset == false || count(ids.begin(),ids.end(),seq.get_id())==1){
                seq.perm_reverse_complement();
            }
            (*poos) << seq.get_fasta();
        }
    }else{
       bool first = true;
       vector<Sequence> done; //for pguess
       while (read_next_seq_from_stream(*pios,ft,retstring,seq)) {
           if (first == true){
               done.push_back(seq);
               (*poos) << seq.get_fasta();
               first = false;
           }else{
               if (reverse_it_or_not(done, seq)){
                    seq.perm_reverse_complement();
               }
               (*poos) << seq.get_fasta();
               if(pguess){
                   done.push_back(seq);
               }else if(sguess){
                   double r = ((double) rand() / (RAND_MAX));
                    if (r < sguess_samplenum){
                        done.push_back(seq);
                    }
               }
           }
        }
        if (ft == 2) {
           if (reverse_it_or_not(done, seq)){
                seq.perm_reverse_complement();
           }
           (*poos) << seq.get_fasta();

        }
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
