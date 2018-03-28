
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstring>
#include <getopt.h>
#include <sstream>

using namespace std;

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "seq_info.h"
#include "log.h"

void print_help() {
    cout << "Print sequence file summary" << endl;
    cout << "By default returns all properties. Alternatively choose 1 property." << endl;
    cout << "This will take fasta, phylip or nexus file formats" << endl;
    cout << endl;
    cout << "Usage: pxlssq [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     input seq file, stdin otherwise" << endl;
    cout << " -i, --indiv         output stats for individual sequences" << endl;
    cout << " -n, --nseq          return the number of sequences" << endl;
    cout << " -c, --nchar         return the number of characters (only if aligned)" << endl;
    cout << "                        - for unaligned seqs, use with -i flag" << endl;
    cout << " -l, --labels        return all taxon labels (one per line)" << endl;
    cout << " -p, --prot          force interpret as protein (if inference fails)" << endl;
    cout << " -a, --aligned       return whether sequences are aligned (same length)" << endl;
    cout << " -f, --freqs         return character state frequencies" << endl;
    cout << " -m, --missing       return the proportion of missing characters" << endl;
    cout << " -o, --outf=FILE     output stats file, stout otherwise" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}
string versionline("pxlssq 0.1\nCopyright (C) 2016 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown and Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"indiv", no_argument, NULL, 'i'},
    {"nseq", no_argument, NULL, 'n'},
    {"nchar", no_argument, NULL, 'c'},
    {"labels", no_argument, NULL, 'l'},
    {"prot", no_argument, NULL, 'p'},
    {"aligned", no_argument, NULL, 'a'},
    {"freqs", no_argument, NULL, 'f'},
    {"missing", no_argument, NULL, 'm'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]){
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    bool indiv = false;
    bool force_protein = false; // i.e. if inference fails
    bool optionsset = false; // is true, do not return all properties
    bool get_labels = false;
    bool check_aligned = false;
    bool get_nseq = false;
    bool get_nchar = false;
    bool get_freqs = false;
    bool get_missing = false;
    char * outf = NULL;
    char * seqf = NULL;
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:o:inclpafmhV", long_options, &oi);
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
            case 'i':
                indiv = true;
                break;
            case 'n':
                get_nseq = true;
                optionsset = true;
                break;
            case 'c':
                get_nchar = true;
                optionsset = true;
                break;
            case 'l':
                get_labels = true;
                optionsset = true;
                break;
            case 'p':
                force_protein = true;
                break;
            case 'a':
                check_aligned = true;
                optionsset = true;
                break;
            case 'f':
                get_freqs = true;
                optionsset = true;
                break;
            case 'm':
                get_missing = true;
                optionsset = true;
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
    
    ostream * poos = NULL;
    ofstream * ofstr = NULL;
    ifstream * fstr = NULL;
    istream * pios = NULL;
    
    if ((get_labels + check_aligned + get_nseq + get_freqs + get_nchar + get_missing) > 1) {
        cout << "Specify 1 property only (or leave blank to show all properties)" << endl;
        exit(0);
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
    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    
    SeqInfo ls_Seq(pios, poos, indiv, force_protein);
    if (optionsset) {
        // get single property
        ls_Seq.get_property (get_labels, check_aligned, get_nseq,
            get_freqs, get_nchar, get_missing);
    } else {
        // the original behaviour
        ls_Seq.summarize();
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

