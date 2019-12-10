#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <algorithm>

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "log.h"
#include "constants.h"

extern std::string PHYX_CITATION;


void print_help() {
    std::cout << "Remove sequences by label" << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus inputs." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxrms [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -s, --seqf=FILE     input nucleotide sequence file, stdin otherwise" << std::endl;
    std::cout << " -n, --names=CSL     names sep by commas (NO SPACES!)" << std::endl;
    std::cout << " -f, --namesf=FILE   names in a file (each on a line)" << std::endl;
    std::cout << " -c, --comp          take the complement (i.e. remove any taxa not in list)" << std::endl;
    std::cout << " -o, --outf=FILE     output sequence file, stout otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxrms 1.0\nCopyright (C) 2015-2020 FePhyFoFum\nLicense GPLv3\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"seqf", required_argument, NULL, 's'},
    {"names", required_argument, NULL, 'n'},
    {"namesf", required_argument, NULL, 'f'},
    {"comp", no_argument, NULL, 'c'},
    {"outf", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
       
    bool fileset = false;
    bool outfileset = false;
    bool namesset = false;
    bool namefileset = false;
    bool complement = false;
    char * namesc = NULL;
    char * namesfc = NULL;
    char * seqf = NULL;
    char * outf = NULL;
    std::string rmf = "";
    std::vector<std::string> names;

    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "s:n:f:co:hVC", long_options, &oi);
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
                std::cout << versionline << std::endl;
                exit(0);
            case 'C':
                std::cout << PHYX_CITATION << std::endl;
                exit(0);
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
    }
    
    if (fileset && outfileset) {
        check_inout_streams_identical(seqf, outf);
    }
    
    std::istream * pios = NULL;
    std::ostream * poos = NULL;
    std::ifstream * fstr = NULL;
    std::ofstream * ofstr = NULL;
    
    if (namesset == true) {
        std::vector<std::string> tokens2;
        std::string del2(",");
        tokens2.clear();
        tokenize(namesc, tokens2, del2);
        for (unsigned int j=0; j < tokens2.size(); j++) {
            trim_spaces(tokens2[j]);
            names.push_back(tokens2[j]);
        }
    } else if (namefileset == true) {
        std::ifstream nfstr(namesfc);
        std::string tline;
        while (getline(nfstr, tline)) {
            trim_spaces(tline);
            names.push_back(tline);
        }
        nfstr.close();
    } else {
        std::cerr << "Error: you need to set the names of the taxa you want to remove (-n). Exiting." << std::endl;
        exit(0);
    }
    
    if (fileset == true) {
        fstr = new std::ifstream(seqf);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (check_for_input_to_stream() == false) {
            print_help();
            exit(1);
        }
    }    
    if (outfileset == true) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
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
                it = find(names.begin(), names.end(), seq_name);
                if ( ((!complement) && (it == names.end())) || ((complement) && (it != names.end())) ) {
                    *poos << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                }
            }
        } else {
            std::vector<Sequence> seqs = read_interleaved_nexus(*pios, num_taxa, num_char);
            for (unsigned int i = 0; i < seqs.size(); i++) {
                seq = seqs[i];
                seq_name = seq.get_id();
                it = find(names.begin(), names.end(), seq_name);
                if ( ((!complement) && (it == names.end())) || ((complement) && (it != names.end())) ) {
                    *poos << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
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
            for (unsigned int i = 0; i < seqs.size(); i++) {
                seq = seqs[i];
                seq_name = seq.get_id();
                it = find(names.begin(), names.end(), seq_name);
                if ( ((!complement) && (it == names.end())) || ((complement) && (it != names.end())) ) {
                    *poos << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                }
            }
        } else {
            // fasta, fastq, or simple phylip
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                seq_name = seq.get_id();
                it = find(names.begin(), names.end(), seq_name);
                if ( ((!complement) && (it == names.end())) || ((complement) && (it != names.end())) ) {
                    *poos << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
                }
            }
            // fasta has a trailing one
            if (ft == 2) {
                seq_name = seq.get_id();
                it = find(names.begin(), names.end(), seq_name);
                if ( ((!complement) && (it == names.end())) || ((complement) && (it != names.end())) ) {
                    *poos << ">" << seq_name << "\n" << seq.get_sequence() << std::endl;
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
