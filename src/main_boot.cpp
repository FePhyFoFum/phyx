/*
 Bare-bones sequence alignment resampling. Default is bootstrap, alternative is jackknife.
 Conserved-partition bootstrap now implemented.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "seq_sample.h"
#include "log.h"

void print_help() {
    cout << "Sequence alignment bootstrap or jackknife resampling." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << "Can read from stdin or file." << endl;
    cout << endl;
    cout << "Usage: pxboot [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input sequence file, stdin otherwise" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise" << endl;
    cout << " -p, --partf=FILE    file listing empirical partitions: NAME = START-STOP[\\INTERVAL]" << endl;
    cout << " -f, --frac=DOUBLE   jackknife percentage, default bootstrap (i.e. -f 1.0)" << endl;
    cout << " -x, --seed=INT      random number seed, clock otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxboot 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv3\nwritten by Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"partf", required_argument, NULL, 'p'},
    {"frac", required_argument, NULL, 'f'},
    {"seed", required_argument, NULL, 'x'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    bool partitioned = false;
    float jackfract = 0.0;
    int numchar = 0;
    char * outf = NULL;
    char * seqf = NULL;
    string partf = "";
    int seed = -1;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:p:f:x:hV", long_options, &oi);
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
            case 'p':
                partitioned = true;
                partf = strdup(optarg);
                break;
            case 'f':
                jackfract = atof(strdup(optarg));
                if (jackfract < 0 || jackfract > 1) {
                    cout << "Jackknife fraction must be 0 < x < 1" << endl;
                    exit(0);
                }
                break;
            case 'x':
                seed = atoi(strdup(optarg));
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
    
    if (partitioned && (jackfract != 0.0)) {
        cout << "Partitioned jackknife not implemented. Exiting." << endl;
        exit (0);
    }
    
    SequenceSampler ss(seed, jackfract, partf);
    
    Sequence seq;
    string retstring;
    bool first = true;
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        if (first) { // need to read in first sequence to get numchar
            numchar = (int)seq.get_sequence().size(); // check against partition information
            if (partitioned) {
                if (numchar != ss.get_num_partitioned_sites()) {
                    cout << "Error: numSites in sequence (" << numchar <<
                        ") does not match that in partition file (" << ss.get_num_partitioned_sites() <<
                        ")." << endl;
                }
            }
            ss.sample_sites(numchar);
            first = false;
        }
        (*poos) << ">" << seq.get_id() << endl;
        (*poos) << ss.get_resampled_seq(seq.get_sequence()) << endl;
    }
    // have to deal with last sequence outside while loop. fix this.
    if (ft == 2) {
        (*poos) << ">" << seq.get_id() << endl;
        (*poos) << ss.get_resampled_seq(seq.get_sequence()) << endl;
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
