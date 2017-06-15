/*
 * main_SW.cpp
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cstring>
#include <getopt.h>

#ifdef OMP
#include <omp.h>
#endif

using namespace std;

#include "utils.h"
#include "seq_reader.h"
#include "sequence.h"
#include "seq_utils.h"
#include "seq_models.h"
#include "pairwise_alignment.h"
#include "log.h"

void print_help() {
    cout << "Conduct Smith-Waterman analysis for all the seqs in a file." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << "Output is a list of the scores and distances (and the alignments" << endl;
    cout << "if asked)." << endl;
    cout << "Can read from stdin or file." << endl;
    cout << endl;
    cout << "Usage: pxsw [OPTION]... [FILE]..." << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output score/distance file, stout otherwise" << endl;
    cout << " -a, --outalnf=FILE  output sequence file, won't output otherwise" << endl;
    cout << " -t, --seqtype=INT   sequence type, default=DNA (DNA=0,AA=1)" << endl;
    cout << " -m, --matrix=FILE   scoring matrix, default DNA=EDNAFULL, AA=BLOSUM62" << endl;
    cout << " -n, --nthreads=INT  number of threads (open mp), default=2" << endl;
    cout << " -v, --verbose       make the output more verbose, turns off parallel" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" <<endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxsw 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"outalnf", required_argument, NULL, 'a'},
    {"seqtype", required_argument, NULL, 't'},
    {"matrix", required_argument, NULL, 'm'},
    {"nthreads", required_argument, NULL, 'n'},
    {"verbose", no_argument, NULL, 'v'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool fileset = false;
    bool outfileset = false;
    bool outalnfileset = false;
    bool matrixfileset = false;
    char * seqf = NULL;
    char * outf = NULL;
    char * outaf = NULL;
    char * matf = NULL;
    int seqtype = 0; //DNA default, 1 = aa
    int num_threads = 2; //DNA default, 1 = aa
    bool verbose = false;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:a:t:m:n:vhV", long_options, &oi);
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
            case 'a':
                outalnfileset = true;
                outaf = strdup(optarg);
                break;
            case 't':
                seqtype = atoi(strdup(optarg));
                if (seqtype > 1) {
                    cout << "Don't recognize seqtype " << seqtype << ". Must be 0 (DNA) or 1 (AA)." << endl;
                    exit(0);
                }
                break;
            case 'm':
                matrixfileset = true;
                matf = strdup(optarg);
                check_file_exists(matf);
                break;
            case 'n':
                num_threads = atoi(strdup(optarg));
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
    map<char, map<char,int> > sc_mat;
    if (matrixfileset == true) {
        read_scoring_matrix(matf, sc_mat);
    } else {
        if (seqtype == 0) {
            get_ednafull(sc_mat);
        } else { //aa
            
            // *** what is supposed to go here? ***
             
        }
    }
    vector<Sequence> seqs;
    Sequence seq;
    string retstring;
    
    istream * pios = NULL;
    ostream * poos = NULL;
    ifstream * fstr = NULL;
    ofstream * ofstr = NULL;
    ofstream * afstr = NULL;
    
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
    if (outalnfileset) {
        afstr = new ofstream(outaf);
    }

    int ft = test_seq_filetype_stream(*pios,retstring);
    while (read_next_seq_from_stream(*pios,ft,retstring,seq)) {
        seqs.push_back(seq);
    }
    // fasta has a trailing one
    if (ft == 2) {
        seqs.push_back(seq);
    }

    // go all by all
    for (unsigned int i=0; i < seqs.size(); i++) {
#ifdef OMP
        omp_set_num_threads(num_threads);
#endif
        #pragma omp parallel for
        for (unsigned int j=0; j < seqs.size(); j++) {
            if (j > i) {
                string aln1;
                string aln2;
                double sc = sw(seqs[i],seqs[j],sc_mat,0, aln1, aln2);
                #pragma omp critical
                {
                    (*poos) << seqs[i].get_id() << "\t" << seqs[j].get_id()  << "\t" << sc << endl;
                    if (verbose) {
                        cout << seqs[i].get_id() <<  "\t" << aln1 << "\n" << seqs[j].get_id()  << "\t" << aln2 << endl;
                    }
                    if (outalnfileset) {
                        (*afstr) << seqs[i].get_id() <<  "\t" << aln1 << "\n" << seqs[j].get_id()  << "\t" << aln2 << endl;
                    }
                }
            }
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
    if (outalnfileset) {
        afstr->close();
    }
    return EXIT_SUCCESS;
}

