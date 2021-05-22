#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstring>
#include <getopt.h>

#include "string_node_object.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "tree_reader.h"
#include "tree.h"
#include "tree_utils.h"
#include "cont_models.h"
#include "optimize_cont_models_nlopt.h"
#include "log.h"
#include "citations.h" // contains PHYX_CITATION


void print_help (void);

void print_help () {
    std::cout << "Continuous character rate estimation with Brownian and OU." << std::endl;
    std::cout << "This will take fasta, phylip, and nexus formats from a file or STDIN." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxcontrates [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -c, --charf=FILE     input character file, STDIN otherwise" << std::endl;
    std::cout << " -t, --treef=FILE     input tree file, STDIN otherwise" << std::endl;
    std::cout << " -a, --analysis=NUM   analysis type (0=anc[DEFAULT], 1=ratetest)" << std::endl;
    std::cout << " -o, --outf=FILE      output sequence file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help           display this help and exit" << std::endl;
    std::cout << " -V, --version        display version and exit" << std::endl;
    std::cout << " -C, --citation       display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

static std::string versionline("pxcontrates 1.2\nCopyright (C) 2013-2021 FePhyFoFum\nLicense GPLv3\nWritten by Stephen A. Smith (blackrim)");

static struct option const long_options[] =
{
    {"char", required_argument, nullptr, 'c'},
    {"tree", required_argument, nullptr, 't'},
    {"outf", required_argument, nullptr, 'o'},
    {"analysis", required_argument, nullptr, 'a'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool cfileset = false;
    bool tfileset = false;
    bool outfileset = false;
    
    char * treef = nullptr;
    char * charf = nullptr;
    char * outf = nullptr;
    int analysis = 0;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "c:t:o:a:hVC", long_options, &oi);
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
            case 'a':
                if (optarg[0] == '1') {
                    analysis = 1;
                }
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
                print_error(argv[0]);
                exit(0);
        }
    }

    std::istream * pios = nullptr;
    std::istream * poos = nullptr;
    std::ifstream * cfstr = nullptr;
    std::ifstream * tfstr = nullptr;

    std::ostream * poouts = nullptr;
    std::ofstream * ofstr = nullptr;
    
    if (tfileset) {
        tfstr = new std::ifstream(treef);
        poos = tfstr;
    } else {
        poos = &std::cin;
    }

    if (cfileset) {
        cfstr = new std::ifstream(charf);
        pios = cfstr;
    } else {
        std::cerr << "Error: you have to set a character file. Only a tree file can be read in through the stream. Exiting." << std::endl;
        exit(1);
    }

    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poouts = ofstr;
    } else {
        poouts = &std::cout;
    }

    std::string retstring;
    int ft = test_char_filetype_stream(*pios, retstring);
    if (ft != 1 && ft != 2) {
        std::cerr << "Error: only fasta and phylip (with spaces) supported so far. Exiting." << std::endl;
        exit(0);
    }
    Sequence seq;
    std::vector<Sequence> seqs;
    std::map<std::string, int> seq_map;
    int y = 0;
    int num_chars = 0;
    while (read_next_seq_char_from_stream(*pios, ft, retstring, seq)) {
        seqs.push_back(seq);
        num_chars = seq.get_num_cont_char();
        seq_map[seq.get_id()] = y;
        seq.clear_cont_char();
        y++;
    }
    if (ft == 2) {
        seqs.push_back(seq);
        seq_map[seq.get_id()] = y;
        seq.clear_cont_char();
    }
    
    // read trees
    TreeReader tr;
    std::vector<Tree *> trees;
    while (getline_safe(*poos, retstring)) {
        if (retstring.empty()) {
            continue;
        }
        trees.push_back(tr.readTree(retstring));
    }
    
    // conduct analyses for each character
    for (int c=0; c < num_chars; c++) {
        std::cerr << "character: " << c << std::endl;
        if (analysis == 0) {
           // std::cout << "Input tree: " << getNewickString(trees[0]) << ";" << std::endl;
            if (c == 0) {
                (*poouts) << "#nexus" << std::endl << "begin trees;" << std::endl;
            }
            for (unsigned int x = 0; x < trees.size(); x++) {
                for (int i = 0; i < trees[x]->getExternalNodeCount(); i++) {
                    std::vector<Superdouble> tv (1);
                    tv[0] = seqs[seq_map[trees[x]->getExternalNode(i)->getName()]].get_cont_char(c);
                    trees[x]->getExternalNode(i)->assocDoubleVector("val", tv);
                }
                for (int i = 0; i < trees[x]->getInternalNodeCount(); i++) {
                    std::vector<Superdouble> tv (1);
                    tv[0] = 0;
                    trees[x]->getInternalNode(i)->assocDoubleVector("val", tv);
                }
                calc_square_change_anc_states(trees[x], 0); // second character dies here
                for (int i = 0; i < trees[x]->getInternalNodeCount(); i++) {
                    double tv = (*trees[x]->getInternalNode(i)->getDoubleVector("val"))[0];
                    trees[x]->getInternalNode(i)->deleteDoubleVector("val");
                    std::ostringstream s;
                    s.precision(9);
                    s << std::fixed << tv;
                    StringNodeObject nob(s.str());
                    trees[x]->getInternalNode(i)->assocObject("value", nob);
                    //trees[x]->getInternalNode(i)->setName(s.str());
                }
                for (int i = 0; i < trees[x]->getExternalNodeCount(); i++) {
                    double tv = (*trees[x]->getExternalNode(i)->getDoubleVector("val"))[0];
                    trees[x]->getExternalNode(i)->deleteDoubleVector("val");
                    std::ostringstream s;
                    s.precision(9);
                    s << std::fixed << tv;
                    StringNodeObject nob(s.str());
                    trees[x]->getExternalNode(i)->assocObject("value", nob);
                    //s << fixed << trees[x]->getExternalNode(i)->getName() << "[&value=" << tv << "]";
                    //trees[x]->getExternalNode(i)->setName(s.str());
                }
                (*poouts) << "tree tree" << c << " = ";
                (*poouts) << getNewickString(trees[x], static_cast<std::string>("value")) << std::endl;
            }
            if (c == (num_chars - 1)) {
                (*poouts) << "end;\n" << std::endl;
            }
            // remove annotations
            remove_annotations(trees[0]);
            
        } else if (analysis == 1) {
            mat vcv;
            int t_ind = 0; // TODO: do this over trees
            int c_ind = c;
            calc_vcv(trees[t_ind], vcv);
            int n = trees[t_ind]->getExternalNodeCount();
            rowvec x = rowvec(n);
            for (int i = 0; i < n; i++) {
                x(i) = seqs[seq_map[trees[t_ind]->getExternalNode(i)->getName()]].get_cont_char(c_ind);
            }
            std::vector<double> res = optimize_single_rate_bm_nlopt(x, vcv, true);
            double aic = (2*2)-(2*(-res[2]));
            double aicc = aic + ((2*2*(2+1))/(n-2-1));
            std::cout << c << " BM " << " state: " << res[0] << " rate: " << res[1]
                << " like: " << -res[2] << " aic: " << aic << " aicc: " << aicc << std::endl;

            std::vector<double> res2 = optimize_single_rate_bm_ou_nlopt(x, vcv);
            aic = (2*3)-(2*(-res2[3]));
            aicc = aic + ((2*3*(3+1))/(n-3-1));
            std::cout << c << " OU " << " state: " << res2[0] << " rate: "
                << res2[1] << " alpha: " << res2[2] << " like: " << -res2[3]
                << " aic: " << aic << " aicc: " << aicc << std::endl;
        }
    }

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
