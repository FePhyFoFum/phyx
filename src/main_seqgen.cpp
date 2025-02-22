#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "seq_gen.h"
#include "utils.h"
#include "tree.h"
#include "tree_reader.h"
#include "log.h"
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Basic sequence simulator under the GTR model." << std::endl;
    std::cout << "This will take fasta, fastq, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxseqgen [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -t, --treef=FILE       input treefile, STDIN otherwise" << std::endl;
    std::cout << " -l, --length=INT       length of sequences to generate. default is 1000" << std::endl;
    std::cout << " -b, --basef=Input      comma-delimited base freqs in order: A,C,G,T. default is equal" << std::endl;
    std::cout << " -g, --gamma=INT        gamma shape value. default is no rate variation" << std::endl;
    std::cout << " -i, --pinvar=FLOAT     proportion of invariable sites. default is 0.0" << std::endl;
    std::cout << " -r, --ratemat=Input    comma-delimited input values for rate matrix. default is JC69" << std::endl;
    std::cout << "                          order: A<->C,A<->G,A<->T,C<->G,C<->T,G<->T" << std::endl;
    std::cout << " -w, --aaratemat=Input  comma-delimited amino acid rate matrix. default is all freqs equal" << std::endl;
    std::cout << "                        order is ARNDCQEGHILKMFPSTWYV" << std::endl;
    std::cout << " -q, --aabasefreq=Input AA frequencies, order: ARNDCQEGHILKMFPSTWYV" << std::endl;
    std::cout << " -c, --protein          run as amino acid" << std::endl;
    std::cout << " -n, --nreps=INT        number of replicates" << std::endl;
    std::cout << " -x, --seed=INT         random number seed, clock otherwise" << std::endl;
    std::cout << " -a, --ancestors        print the ancestral node sequences. default is no" << std::endl;
    std::cout << "                          use -p for the nodes labels" << std::endl;
    std::cout << " -p, --printnodelabels  print newick with internal node labels. default is no" << std::endl;
    std::cout << " -m, --multimodel=Input specify multiple models across tree" << std::endl;
    std::cout << "                          input is as follows:" << std::endl;
    std::cout << "                            A<->C,A<->G,A<->T,C<->G,C<->T,G<->T,Node#,A<->C,A<->G,A<->T,C<->G,C<->T,G<->T" << std::endl;
    std::cout << "                            EX:.3,.3,.3,.3,.3,1,.3,.3,.2,.5,.4" << std::endl;
    std::cout << " -k, --rootseq=STRING   set root sequence. default is random (from basefreqs)" << std::endl;
    std::cout << " -o, --outf=FILE        output seq file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help             display this help and exit" << std::endl;
    std::cout << " -V, --version          display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string get_version_line () {
    std::string vl = "pxseqgen 1.3.2\n";
    vl += "Copyright (C) 2015-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)";
    return vl;
}

static struct option const long_options[] =
{
    {"treef", required_argument, nullptr, 't'},
    {"outf", required_argument, nullptr, 'o'},
    {"length", required_argument, nullptr, 'l'},
    {"basef", required_argument, nullptr, 'b'},
    {"gamma", required_argument, nullptr, 'g'},
    {"pinvar", required_argument, nullptr, 'i'},
    {"ratemat", required_argument, nullptr, 'r'},
    {"aaratemat", required_argument, nullptr, 'w'},
    {"aabasef", required_argument, nullptr, 'q'},
    {"nreps", required_argument, nullptr, 'n'},
    {"seed", required_argument, nullptr, 'x'},
    {"ancestors", no_argument, nullptr, 'a'},
    {"printnodelabels", no_argument, nullptr, 'p'},
    {"protein", no_argument, nullptr, 'c'},
    {"multimodel", required_argument, nullptr, 'm'},
    {"rootseq", required_argument, nullptr, 'k'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool fileset = false;
    bool printpost = false;
    bool showancs = false;
    bool is_dna = true;
    float pinvar = 0.0;
    double tot;
    std::string yorn = "n";
    unsigned int seqlen = 1000;
    int x = 0;
    int pos = 0;
    int pos2 = 0;
    std::string infreqs;
    std::string inrates;
    std::string holdrates;
    std::string ancseq;
    char * outf = nullptr;
    char * treef = nullptr;
    std::vector<double> diag(20, 0.0);
    std::vector<double> basefreq(4, 0.25);
    std::vector<double> aabasefreq(20, 0.05);
    std::vector<double> userrates;
    std::vector<double> multirates;
    int nreps = 1; // not implemented at the moment
    long int seed = -1;
    int numpars = 0;
    float alpha = -1.0;
    std::vector< std::vector<double>> dmatrix;
    std::vector< std::vector<double> > aa_rmatrix(20, std::vector<double>(20, 1));
        for (unsigned int i = 0; i < aa_rmatrix.size(); i++) {
        for (unsigned int j = 0; j < aa_rmatrix.size(); j++) {
            if (i == j) { // Fill Diagonal
                aa_rmatrix[i][j] = -19.0;
            }
        }
    }
    std::vector< std::vector<double> > rmatrix(4, std::vector<double>(4, 0.33));
    for (unsigned int i = 0; i < rmatrix.size(); i++) {
        for (unsigned int j = 0; j < rmatrix.size(); j++) {
            if (i == j) { // Fill Diagonal
                rmatrix[i][j] = -0.99;
            }
        }
    }
    /*dmatrix = aa_rmatrix;
    for (unsigned int i = 0; i < dmatrix.size(); i++) {
        for (unsigned int j = 0; j < dmatrix.size(); j++) {
            std::cout << dmatrix[i][j] << " ";
        }
        std::cout << "\n";
    vl += "";
    }*/

    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "t:o:l:b:g:i:r:w:q:n:x:apcm:k:hVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 't':
                fileset = true;
                treef = strdup(optarg);
                check_file_exists(treef);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'b':
                infreqs = strdup(optarg);
                parse_comma_list(infreqs, basefreq);
                if (basefreq.size() != 4) {
                    std::cerr << "Error: must provide 4 base frequencies (" << basefreq.size()
                        << " provided). Exiting." << std::endl;
                    exit(0);
                }
                if (!essentially_equal(sum(basefreq), 1.0)) {
                    std::cerr << "Error: base frequencies must sum to 1.0. Exiting." << std::endl;
                    exit(0);
                }
                break;
            case 'l':
                x = string_to_int(optarg, "-l");
                if (x < 0) {
                    std::cerr << "Error: Sequence length must be a positive integer. Exiting."
                            << std::endl;
                    exit(0);
                }
                seqlen = static_cast<unsigned int>(x);
                break;
            case 'a':
                showancs = true;
                break;
            case 'r':
                inrates = strdup(optarg);
                parse_comma_list(inrates, userrates);
                
                // NOTE: will have to alter this check for a.a., non-reversible, etc.
                if (userrates.size() != 6) {
                    std::cerr << "Error: must provide 6 substitution parameters. " <<
                        "Only " << userrates.size() << " provided. Exiting." << std::endl;
                    exit(0);
                }
                // NOTE: this uses order: A,T,C,G for matrix, but
                //       A<->C,A<->G,A<->T,C<->G,C<->T,G<->T for subst. params.
                rmatrix[0][2] = userrates[0];
                rmatrix[2][0] = userrates[0];
                rmatrix[0][3] = userrates[1];
                rmatrix[3][0] = userrates[1];
                rmatrix[0][1] = userrates[2];
                rmatrix[1][0] = userrates[2];
                rmatrix[2][3] = userrates[3];
                rmatrix[3][2] = userrates[3];
                rmatrix[1][2] = userrates[4];
                rmatrix[2][1] = userrates[4];
                rmatrix[1][3] = userrates[5];
                rmatrix[3][1] = userrates[5];
                rmatrix[0][0] = (userrates[0]+userrates[1]+userrates[2]) * -1;
                rmatrix[1][1] = (userrates[2]+userrates[4]+userrates[5]) * -1;
                rmatrix[2][2] = (userrates[0]+userrates[3]+userrates[4]) * -1;
                rmatrix[3][3] = (userrates[1]+userrates[3]+userrates[5]) * -1;
                /*//Turn on to check matrix
                for (unsigned int i = 0; i < rmatrix.size(); i++) {
                   for (unsigned int j = 0; j < rmatrix.size(); j++) {
                      std::cout << rmatrix[i][j] << " ";
                   }
                    std::cout << "\n";
    vl += "";
                }*/
                break;
            case 'w':
                inrates = strdup(optarg);
                parse_comma_list(inrates, userrates);
                is_dna = false;
                
                // NOTE: will have to alter this check for a.a., non-reversible, etc.
                if (userrates.size() != 190) {
                    std::cerr << "Error: must provide 190 substitution parameters. ";
                    std::cerr << "I know its a stupidly large amount. " <<
                        "Only " << userrates.size() << " provided. Exiting." << std::endl;
                    exit(0);
                }
                pos = 0;
                pos2 = 1;
                //Fill the Matrix
                for (unsigned int i = 0; i < userrates.size(); i++) {
                    aa_rmatrix[pos][pos2] = userrates[i];
                    aa_rmatrix[pos2][pos] = userrates[i];
                    pos2++;
                    if (pos2 == 20) {
                        pos += 1;
                        pos2 = (pos + 1);
                    }
                }
                //Replace Diagonal
                for (unsigned int i = 0; i < aa_rmatrix.size(); i++) {
                    for (unsigned int j = 0; j < aa_rmatrix.size(); j++) {
                        if (i != j) {
                            tot += aa_rmatrix[i][j];
                        }
                    }
                    aa_rmatrix[i][i] = (tot*-1);
                    tot = 0.0;
                }
                /*
                for (unsigned int i = 0; i < aa_rmatrix.size(); i++) {
                    for (unsigned int j = 0; j < aa_rmatrix.size(); j++) {
                        std::cout << aa_rmatrix[i][j] << " ";
                    }
                    std::cout << "\n";
    vl += "";
                }*/
                break;
            case 'n':
                nreps = string_to_int(optarg, "-n");
                break;
            case 'x':
                seed = string_to_long_int(optarg, "-x");
                break;
            case 'q':
                is_dna = false;
                infreqs = strdup(optarg);
                parse_comma_list(infreqs, aabasefreq);
                if (aabasefreq.size() != 20) {
                    std::cerr << "Error: must provide 20 base frequencies (" << aabasefreq.size()
                        << " provided). Exiting." << std::endl;
                    exit(0);
                }
                if (!essentially_equal(sum(aabasefreq), 1.0)) {
                    std::cerr << "Error: base frequencies must sum to 1.0. Exiting." << std::endl;
                    exit(0);
                }
                break;
            case 'g':
                alpha = string_to_float(optarg, "-g");
                break;
            case 'i':
                pinvar = string_to_float(optarg, "-i");
                break;
            case 'p':
                printpost = true;
                break;
            case 'c':
                is_dna = false;
                break;
            case 'm':
                holdrates = strdup(optarg);
                parse_comma_list(holdrates, multirates);
                numpars = multirates.size();
                if ((numpars - 6) % 7 != 0) {
                    std::cerr << "Error: must provide 6 background substitution "
                        << "parameters and 7 values (1 node id + 6 subst. par.) "
                        << "for each piecewise model. Exiting." << std::endl;
                    exit(0);
                }
                break;
            case 'k':
                ancseq = strdup(optarg);
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
        check_inout_streams_identical(treef, outf);
    }
    
    if (is_dna) {
        dmatrix = rmatrix;
    } else {
        dmatrix = aa_rmatrix;
    }
    
    std::istream * pios = nullptr;
    std::ostream * poos = nullptr;
    std::ifstream * fstr = nullptr;
    std::ofstream * ofstr = nullptr;
    
    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    if (fileset) {
        fstr = new std::ifstream(treef);
        pios = fstr;
    } else {
        pios = &std::cin;
        if (!check_for_input_to_stream()) {
            print_help();
            exit(1);
        }
    }
    
    
    /*
     * Default Base Frequencies and Rate Matrix
     *
     */
    //vector<double> basefreq(4, 0.0);
    //basefreq[0] = .25;
    //basefreq[1] = .25;
    //basefreq[2] = .25;
    //basefreq[3] = 1.0 - basefreq[0] - basefreq[1] - basefreq[2];
    /*    
    vector< vector<double> > rmatrix(4, vector<double>(4, 0.33));
    for (unsigned int i = 0; i < rmatrix.size(); i++) {
        for (unsigned int j = 0; j < rmatrix.size(); j++) {
            if (i == j) {//Fill Diagnol
                rmatrix[i][j] = -1.0;
            }
        }
    }
    */
    std::string retstring;
    int ft = test_tree_filetype_stream(*pios, retstring);
    if (ft != 0 && ft != 1) {
        std::cerr << "Error: this really only works with nexus or newick. Exiting." << std::endl;
        exit(0);
    }
    
    // allow > 1 tree in input. passing but not yet using nreps
    int treeCounter = 0;
    bool going = true;
    if (ft == 1) { // newick. easy
        while (going) {
            Tree * tree = read_next_tree_from_stream_newick (*pios, retstring, &going);
            if (tree != nullptr) {
                //std::cout << "Working on tree #" << treeCounter << std::endl;
                SequenceGenerator SGen(seqlen, basefreq, dmatrix, tree, showancs,
                    nreps, seed, alpha, pinvar, ancseq, printpost, multirates, aabasefreq, is_dna);
                std::vector<Sequence> seqs = SGen.get_sequences();
                for (unsigned int i = 0; i < seqs.size(); i++) {
                    Sequence seq = seqs[i];
                    (*poos) << ">" << seq.get_id() << std::endl;
                    //std::cout << "Here" << std::endl;
                    (*poos) << seq.get_sequence() << std::endl;
                }
                delete tree;
                treeCounter++;
            }
        }
    } else if (ft == 0) { // Nexus. need to worry about possible translation tables
        std::map<std::string, std::string> translation_table;
        bool ttexists;
        ttexists = get_nexus_translation_table(*pios, &translation_table, &retstring);
        while (going) {
            Tree * tree = read_next_tree_from_stream_nexus(*pios, retstring, ttexists,
                &translation_table, &going);
            if (going) {
                //std::cout << "Working on tree #" << treeCounter << std::endl;
                SequenceGenerator SGen(seqlen, basefreq, dmatrix, tree, showancs,
                    nreps, seed, alpha, pinvar, ancseq, printpost, multirates, aabasefreq, is_dna);
                std::vector<Sequence> seqs = SGen.get_sequences();
                for (unsigned int i = 0; i < seqs.size(); i++) {
                    Sequence seq = seqs[i];
                    (*poos) << ">" << seq.get_id() << std::endl;
                    (*poos) << seq.get_sequence() << std::endl;
                }
                delete tree;
                treeCounter++;
            }
        }
    }
    return EXIT_SUCCESS;
}
