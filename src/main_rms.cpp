#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <algorithm>
#include <regex>

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Remove sequences by label." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << "Results are written in fasta format." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxrms [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s, --seqf=FILE     input nucleotide sequence file, STDIN otherwise" << std::endl;
    std::cout << " -n, --names=CSL     names sep by commas (NO SPACES!)" << std::endl;
    std::cout << " -f, --namesf=FILE   names in a file (each on a line)" << std::endl;
    std::cout << " -r, --regex=STRING  match tip labels by a regular expression" << std::endl;
    std::cout << " -c, --comp          take the complement (i.e. remove any taxa not in list)" << std::endl;
    std::cout << " -o, --outf=FILE     output sequence file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxrms 1.3.1\n";
    vl += "Copyright (C) 2015-2024 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph W. Brown, Joseph F. Walker, Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"seqf", required_argument, nullptr, 's'},
    {"names", required_argument, nullptr, 'n'},
    {"namesf", required_argument, nullptr, 'f'},
    {"regex", required_argument, nullptr, 'r'},
    {"comp", no_argument, nullptr, 'c'},
    {"outf", required_argument, nullptr, 'o'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
       
    bool fileset = false;
    bool outfileset = false;
    bool namesset = false;
    bool namefileset = false;
    bool complement = false;
    bool regex = false;
    std::regex regexp;
    std::string regex_pattern;
    bool match = false; // for regex searches
    
    char * namesc = nullptr;
    char * namesfc = nullptr;
    char * seqf = nullptr;
    char * outf = nullptr;
    std::string rmf;
    std::vector<std::string> names;

    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:n:f:r:co:hVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 's':
                fileset = true;
                seqf = strdup(optarg);
                check_file_exists(seqf);
                break;
            case 'n':
                namesset = true;
                namesc = strdup(optarg);
                break;
            case 'f':
                namefileset = true;
                namesfc = strdup(optarg);
                check_file_exists(namesfc);
                break;
            case 'r':
                regex = true;
                regex_pattern = strdup(optarg);
                regexp.assign(regex_pattern);
                break;
            case 'c':
                complement = true;
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
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
                print_error(*argv);
                exit(0);
        }
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
    if (fileset) {
        fstr = new std::ifstream(seqf);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (!check_for_input_to_stream()) {
            print_help();
            exit(1);
        }
    }    
    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    if (namesset) {
        std::vector<std::string> tokens2;
        std::string del2(",");
        tokens2.clear();
        tokenize(namesc, tokens2, del2);
        for (auto & tk : tokens2) {
            trim_spaces(tk); // this will never have to be used, as spaces would break cmd line call
            names.push_back(tk);
        }
    } else if (namefileset) {
        std::ifstream nfstr(namesfc);
        std::string tline;
        while (getline_safe(nfstr, tline)) {
            trim_spaces(tline);
            if (tline.empty()) {
                continue;
            }
            names.push_back(tline);
        }
        nfstr.close();
    } else if (!regex) {
        std::cerr << "Error: you must specify which tips to remove." << std::endl;
        std::cerr << "This can be done with a list (-n) or file (-f) of names, or a regular expression (-r)."
                << std::endl;
        std::cerr << "Exiting." << std::endl;
        exit(0);
    }
    
    Sequence seq;
    std::string retstring;
    std::string seq_name;
    int num_taxa, num_char; // not used, but required by some readers
    
    int ft = test_seq_filetype_stream(*pios, retstring);
    std::vector<std::string>::iterator it;
    
    // extra stuff to deal with possible interleaved nexus
    if (ft == 0) {
        bool interleave = false;
        get_nexus_dimensions(*pios, num_taxa, num_char, interleave);
        retstring = ""; // need to do this to let seqreader know we are mid-file
        if (!interleave) {
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                seq_name = seq.get_id();
                if (regex) {
                    match = std::regex_search(seq_name, regexp);
                    if ( (match && complement) || (!match && !complement) ) {
                        (*poos) << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                    }
                } else {
                    it = find(names.begin(), names.end(), seq_name);
                    if ( ((!complement) && (it == names.end())) || ((complement) && (it != names.end())) ) {
                        (*poos) << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                    }
                }
            }
        } else {
            std::vector<Sequence> seqs = read_interleaved_nexus(*pios, num_taxa, num_char);
            for (const auto & sq : seqs) {
                seq = sq;
                seq_name = seq.get_id();
                if (regex) {
                    match = std::regex_search(seq_name, regexp);
                    if ( (match && complement) || (!match && !complement) ) {
                        (*poos) << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                    }
                } else {
                    it = find(names.begin(), names.end(), seq_name);
                    if ( ((!complement) && (it == names.end())) || ((complement) && (it != names.end())) ) {
                        (*poos) << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                    }
                }
            }
        }
    } else {
        bool complicated_phylip = false;
        // check if we are dealing with a complicated phylip format
        if (ft == 1) {
            get_phylip_dimensions(retstring, num_taxa, num_char);
            complicated_phylip = is_complicated_phylip(*pios, num_char);
        }
        if (complicated_phylip) {
            std::vector<Sequence> seqs = read_phylip(*pios, num_taxa, num_char);
            for (const auto & sq : seqs) {
                seq = sq;
                seq_name = seq.get_id();
                if (regex) {
                    match = std::regex_search(seq_name, regexp);
                    if ( (match && complement) || (!match && !complement) ) {
                        (*poos) << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                    }
                } else {
                    it = find(names.begin(), names.end(), seq_name);
                    if ( ((!complement) && (it == names.end())) || ((complement) && (it != names.end())) ) {
                        (*poos) << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                    }
                }
            }
        } else {
            // fasta, fastq, or simple phylip
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                seq_name = seq.get_id();
                if (regex) {
                    match = std::regex_search(seq_name, regexp);
                    if ( (match && complement) || (!match && !complement) ) {
                        (*poos) << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                    }
                } else {
                    it = find(names.begin(), names.end(), seq_name);
                    if ( ((!complement) && (it == names.end())) || ((complement) && (it != names.end())) ) {
                        (*poos) << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                    }
                }
            }
            // fasta has a trailing one
            if (ft == 2) {
                seq_name = seq.get_id();
                if (regex) {
                    match = std::regex_search(seq_name, regexp);
                    if ( (match && complement) || (!match && !complement) ) {
                        (*poos) << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                    }
                } else {
                    it = find(names.begin(), names.end(), seq_name);
                    if ( ((!complement) && (it == names.end())) || ((complement) && (it != names.end())) ) {
                        (*poos) << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                    }
                }
            }
        }
    }
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
