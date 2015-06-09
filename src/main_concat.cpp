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
#include <string.h>
#include <getopt.h>

using namespace std;

#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "concat.h"

void print_help() {
    cout << "Sequence file concatenation." << endl;
    cout << "This will take fasta, fastq, phylip, and nexus inputs." << endl;
    cout << "Individual files can be of different formats." << endl;
    cout << endl;
    cout << "Usage: pxconcat [OPTION]... " << endl;
    cout << endl;
    cout << " -s, --seqf=FILE     list of input sequence files" << endl;
    cout << " -o, --outf=FILE     output sequence file, stout otherwise" << endl;
    cout << " -p, --partf=FILE    output partition file, none otherwise" << endl;
    cout << "     --help          display this help and exit" << endl;
    cout << "     --version       display version and exit" << endl;
    cout << endl;
    cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << endl;
    cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << endl;
}

string versionline("pxconcat 0.1\nCopyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"outf", required_argument, NULL, 'o'},
    {"partf", required_argument, NULL, 'p'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    bool outfileset = false;
    bool fileset = false;
    bool logparts = false;
    vector <string> inputFiles;
    SequenceConcatenater result;
    char * outf;
    char * seqf;
    string partf = "";
    int seed = -1;

    while (1) {
        int oi = -1;
        int curind = optind;
        int c = getopt_long(argc, argv, "s:o:p:hV", long_options, &oi);
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
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'p':
                logparts = true;
                partf = strdup(optarg);
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

    if (!fileset || inputFiles.size() < 2) {
        cout << "Must specify 2 or more files to concatenate. Exiting." << endl;
        exit(0);
    }

    ostream* poos;
    ofstream* ofstr;

    if (outfileset == true) {
        ofstr = new ofstream(outf);
        poos = ofstr;
    } else {
        poos = &cout;
    }

    bool first = true;

    for (int i = 0; i < (int)inputFiles.size(); i++) {
        SequenceConcatenater curr(inputFiles[i]);
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
    result.write_partition_information(partf);
    }

    return EXIT_SUCCESS;
}

