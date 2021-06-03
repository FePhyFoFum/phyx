#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "concat.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << "Sequence file concatenation." << std::endl;
    std::cout << "Can use wildcards e.g.:" << std::endl;
    std::cout << "  pxcat -s *.phy -o my_cat_file.fa" << std::endl;
    std::cout << "However, if the argument list is too long (shell limit), put filenames in a file:" << std::endl;
    std::cout << "  for x in *.phy" << std::endl;
    std::cout << "    do echo $x >> flist.txt" << std::endl;
    std::cout << "  done" << std::endl;
    std::cout << "and call using the -f option:" << std::endl;
    std::cout << "  pxcat -f flist.txt -o my_cat_file.fa" << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus sequence formats." << std::endl;
    std::cout << "Individual files may be of different formats." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxcat [OPTIONS]... FILES" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     list of input sequence files (space delimited)" << std::endl;
    std::cout << " -f, --flistFILE     file listing input files (one per line)" << std::endl;
    std::cout << " -p, --partf=FILE    output partition file, none otherwise" << std::endl;
    std::cout << " -u, --uppercase     export characters in uppercase" << std::endl;
    std::cout << " -o, --outf=FILE     output sequence file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxcat 1.2\n";
    vl += "Copyright (C) 2015-2021 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"flist", required_argument, nullptr, 'f'},
    {"outf", required_argument, nullptr, 'o'},
    {"partf", required_argument, nullptr, 'p'},
    {"uppercase", no_argument, nullptr, 'u'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    bool logparts = false;
    bool toupcase = false;
    std::vector<std::string> inputFiles;
    char * outf = nullptr;
    std::string partf;
    std::string listf;

    while (true) {
        int oi = -1;
        int curind;
        int c = getopt_long(argc, argv, "s:f:o:p:uhVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                curind = optind - 1;
                while (curind < argc) {
                    std::string temp = strdup(argv[curind]);
                    curind++;
                    if (temp[0] != '-') {
                        std::ifstream infile(temp.c_str());
                        if (infile.good()) { // check that file exists
                            inputFiles.push_back(temp);
                            infile.close();
                        } else {
                            std::cerr << "Error: cannot find input file '" << temp << "'. Exiting." << std::endl;
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
                std::cout << get_version_line() << std::endl;
                exit(0);
            case 'C':
                std::cout << get_phyx_citation() << std::endl;
                exit(0);
            default:
                print_error(argv[0]);
                exit(0);
        }
    }
    
    if (!fileset) {
        std::cerr << "Error: must specify 1 or more files to concatenate. Exiting." << std::endl;
        exit(0);
    }
    if (!listf.empty()) {
        std::string line;
        std::ifstream ifs(listf.c_str());
        while (getline (ifs, line)) {
            if (!line.empty()) {
                inputFiles.push_back(line);
            }
        }
        ifs.close();
    }
    
    std::ostream * poos = nullptr;
    std::ofstream * ofstr = nullptr;

    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    SequenceConcatenater result(toupcase);
    bool first = true;

    for (auto & inputFile : inputFiles) {
        SequenceConcatenater curr(inputFile, toupcase);
        if (!first) {
            result.concatenate(curr);
        } else {
            result = curr;
            first = false;
        }
    }

    // write sequences
    for (unsigned int i = 0; i < result.get_num_taxa(); i++) {
        Sequence curr = result.get_sequence(i);
        (*poos) << ">" << curr.get_id() << std::endl;
        (*poos) << curr.get_sequence() << std::endl;
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
