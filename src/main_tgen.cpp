#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

#include "tgen.h"
#include "utils.h"
#include "tree_utils.h"
#include "log.h"
#include "constants.h"

extern std::string PHYX_CITATION;


void print_help() {
    std::cout << "Generate all tree topologies for n taxa." << std::endl;
    std::cout << "Random tree samples are a-coming." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxtgen [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -n, --ntax=INT      number of taxa" << std::endl;
    std::cout << " -r, --rooted        whether generated trees are rooted (default: false)" << std::endl;
    std::cout << " -c, --count         give the number of possible trees for n taxa and exit" << std::endl;
    std::cout << " -l, --label=STRING  prefix label for taxon names (default: 't')" << std::endl;
    std::cout << " -o, --outf=FILE     output file, STOUT otherwise" << std::endl;
//    std::cout << " -x, --seed=INT      random number seed, clock otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

std::string versionline("pxtgen 1.0\nCopyright (C) 2020 FePhyFoFum\nLicense GPLv3\nWritten by Joseph W. Brown");

static struct option const long_options[] =
{
    {"ntax", required_argument, NULL, 'n'},
    {"rooted", no_argument, NULL, 'r'},
    {"count", no_argument, NULL, 'r'},
    {"label", required_argument, NULL, 'l'},
    {"outf", required_argument, NULL, 'o'},
//    {"seed", required_argument, NULL, 'x'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {"citation", no_argument, NULL, 'C'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    int num_taxa = 0;
    bool rooted = false;
    bool count = false;
    bool outfileset = false;
    std::string lprefix = "t";
    char * outf = NULL;
    
    // limit on number of terminals supported (exhaustive)
    int sim_limit_exh = 10;
    // bc of the way trees are simulated, cannot be arbitrarily large
    // hope to fix this soon
    
//    int seed = -1;
    
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "n:rcl:o:hVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'n':
                num_taxa = string_to_int(optarg, "-n");
                break;
            case 'r':
                rooted = true;
                break;
            case 'c':
                count = true;
                break;
            case 'l':
                lprefix = strdup(optarg);
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
//            case 'x':
//                seed = string_to_int(optarg, "-x");
//                break;
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
    
    std::string rootstat = (rooted) ? "rooted" : "unrooted";
    
    if (num_taxa == 0) {
        std::cerr << "Error: you have to set the number of taxa -n. Exiting." << std::endl;
        exit(0);
    } else if (num_taxa < 3) {
        std::cerr << "Error: the number of taxa -n must be >= 3. Exiting." << std::endl;
        exit(0);
    } else if (num_taxa > sim_limit_exh) {
        std::cerr << "Error: the number of taxa -n is currently limited to " << sim_limit_exh
                << " (" << get_num_possible_trees(sim_limit_exh, rooted) << " "
                << rootstat << " topologies). Exiting." << std::endl;
        exit(0);
    }
    
    if (count) {
        std::cout << "There are " << get_num_possible_trees(num_taxa, rooted)
                << " possible " << rootstat << " topologies for " << num_taxa
                << " taxa." << std::endl;
        exit(0);
    }
    
    std::ostream * poos = NULL;
    std::ofstream * ofstr = NULL;
    
    if (outfileset == true) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    TopologyGenerator TG(num_taxa, rooted, lprefix);
    TG.get_newicks(poos);
    
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
