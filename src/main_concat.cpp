/*
 * Main_concat.cpp
 *
 *  Created on: Sep 22, 2014
 *      Author: joe
*/

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "concat.h"
#include "log.h"

void print_help() {
    cout << "Sequence file concatenation." << endl;
    cout << "Can use wildcards e.g.:" << endl;
    cout << "  pxcat -s *.phy -o my_cat_file.fa" << endl;
    cout << "However, if the argument list is too long (shell limit), put filenames in a file:" << endl;
    cout << "  for x in *.phy; do echo $x >> flist.txt; done" << endl;
    cout << "and call using the -f option:" << endl;
    cout << "  pxcat -f flist.txt -o my_cat_file.fa" << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << "Individual files can be of different formats." << endl;
    cout << endl;
    cout << "Usage: pxcat [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     list of input sequence files (space delimited)" << endl;
    cout << " -f, --flistFILE     file listing input files (one per line)" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise" << endl;
    cout << " -p, --partf=FILE    output partition file, none otherwise" << endl;
    cout << " -u, --uppercase     export characters in uppercase" << endl;
    cout << " -h, --help          display this help and exit" << endl;
    cout << " -V, --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxcat 0.9\nCopyright (C) 2019 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"flist", required_argument, NULL, 'f'},
    {"outf", required_argument, NULL, 'o'},
    {"partf", required_argument, NULL, 'p'},
    {"uppercase", no_argument, NULL, 'u'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    bool logparts = false;
    bool toupcase = false;
    vector <string> inputFiles;
    char * outf = NULL;
    string partf = "";
    string listf = "";

    while (1) {
        int oi = -1;
        int curind = optind;
        int c = getopt_long(argc, argv, "s:f:o:p:uhV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                curind = optind - 1;
                while (curind < argc) {
                    string temp = strdup(argv[curind]);
                    curind++;
                    if (temp[0] != '-') {
                        ifstream infile(temp.c_str());
                        if (infile.good()) { // check that file exists
                            inputFiles.push_back(temp);
                            infile.close();
                        } else {
                            cout << "Cannot find input file '" << temp << "'. Exiting." << endl;
                            exit(0);
                        }
                    } else {
                        optind = curind - 1;
                        break;
                    }
                }
                break;
            case 'f':
                fileset = true;
                listf = strdup(optarg);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'p':
                logparts = true;
                partf = strdup(optarg);
                break;
            case 'u':
                toupcase = true;
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
    
    if (!fileset) {
        cout << "Must specify 1 or more files to concatenate. Exiting." << endl;
        exit(0);
    }
    if (listf != "") {
        string line;
        ifstream ifs(listf.c_str());
        while (getline (ifs, line)) {
            if (!line.empty()) {
                inputFiles.push_back(line);
            }
        }
        ifs.close();
    }
    
    ostream * poos = NULL;
    ofstream * ofstr = NULL;

    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }
    
    SequenceConcatenater result;
    bool first = true;

    for (unsigned int i = 0; i < inputFiles.size(); i++) {
        SequenceConcatenater curr(inputFiles[i], toupcase);
        if (!first) {
            result.concatenate(curr);
        } else {
            result = curr;
            first = false;
        }
    }

    // write sequences. currently only fasta format.
    for (int i = 0; i < result.get_num_taxa(); i++) {
        Sequence curr = result.get_sequence(i);
        (*poos) << ">" << curr.get_id() << endl;
        (*poos) << curr.get_sequence() << endl;
    }

    if (outfileset) {
        ofstr->close();
        delete poos;
    }

    // log partition information. currently only RAxML-style.
    if (logparts) {
        result.write_partition_information(inputFiles, partf);
    }

    return EXIT_SUCCESS;
}

