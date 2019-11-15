#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>
#include <map>

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "tree_reader.h"
#include "tree.h"
#include "cont_models.h"
#include "optimize_cont_models_nlopt.h"
#include "log.h"

void print_help() {
    std::cout << "Continuous character branch length estimation with Brownian." << std::endl;
    std::cout << "This will take fasta, phylip (and soon nexus) inputs." << std::endl;
    std::cout << "Can read from stdin or file." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxcontrates [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -c, --charf=FILE     input character file, stdin otherwise" << std::endl;
    std::cout << " -t, --treef=FILE     input tree file, stdin otherwise" << std::endl;
    std::cout << " -o, --outf=FILE      output sequence file, stout otherwise" << std::endl;
    std::cout << " -h, --help           display this help and exit" << std::endl;
    std::cout << " -V, --version        display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

string versionline("pxcontbl 0.1\nCopyright (C) 2017 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"char", required_argument, NULL, 'c'},
    {"tree", required_argument, NULL, 't'},
    {"outf", required_argument, NULL, 'o'},
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
        int c = getopt_long(argc, argv, "c:t:o:hV", long_options, &oi);
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
            case 'h':
                print_help();
                exit(0);
            case 'V':
                std::cout << versionline << std::endl;
                exit(0);
            default:
                print_error(argv[0], (char)c);
                exit(0);
        }
    }

    std::istream * pios = NULL;
    std::istream * poos = NULL;
    std::ifstream * cfstr = NULL;
    std::ifstream * tfstr = NULL;

    std::ostream * poouts = NULL;
    std::ofstream * ofstr = NULL;
    

    if (tfileset == true) {
        tfstr = new ifstream(treef);
        poos = tfstr;
    } else {
        poos = &std::cin;
    }

    if (cfileset == true) {
        cfstr = new std::ifstream(charf);
        pios = cfstr;
    } else {
        std::cout << "you have to set a character file. Only a tree file can be read in through the stream;" << std::endl;
    }

    //out file
    //
    if (outfileset == true){
        ofstr = new std::ofstream(outf);
        poouts = ofstr;
    } else{
        poouts = &std::cout;
    }
    //

    std::string retstring;
    int ft = test_char_filetype_stream(*pios, retstring);
    if (ft != 1 && ft != 2) {
        std::cout << "only fasta and phylip (with spaces) supported so far" << std::endl;
        exit(0);
    }
    Sequence seq;
    std::vector<Sequence> seqs;
    std::map<string, int> seq_map;
    int y = 0;
    int nchars = 0 ;
    while (read_next_seq_char_from_stream(*pios, ft, retstring, seq)) {
        seqs.push_back(seq);
        nchars = seq.get_num_cont_char();
        seq_map[seq.get_id()] = y;
        seq.clear_cont_char();
        y++;
    }
    std::cout << "nchars: " <<  nchars << std::endl;
    
    if (ft == 2) {
        seqs.push_back(seq);
        seq_map[seq.get_id()] = y;
        seq.clear_cont_char();
    }
    //read trees
    TreeReader tr;
    std::vector<Tree *> trees;
    while (getline(*poos,retstring)) {
        if (retstring.size()<4){
            continue;
        }
        trees.push_back(tr.readTree(retstring));
    }
    int x = 0;
    //conduct analyses for each character
    for (int i=0; i < trees[x]->getExternalNodeCount(); i++) {
        std::vector<Superdouble> tv (nchars);
        for (int c=0; c < nchars; c++) {
            tv[c] = seqs[seq_map[trees[x]->getExternalNode(i)->getName()]].get_cont_char(c);
        }
        trees[x]->getExternalNode(i)->assocDoubleVector("val",tv);
    }
    for (int i=0; i < trees[x]->getInternalNodeCount(); i++) {
        std::vector<Superdouble> tv (nchars);
        for (int c=0; c < nchars; c++) {
            tv[c] = 0;
        }
        trees[x]->getInternalNode(i)->assocDoubleVector("val",tv);
    }
    float sigma = 1;
    std::cout << calc_bm_prune(trees[x], sigma) << std::endl;
    optimize_single_rate_bm_bl(trees[x]);
    (*poouts) << trees[x]->getRoot()->getNewick(true) << ";" << std::endl;
    std::cout << calc_bm_prune(trees[x], sigma) << std::endl;

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
