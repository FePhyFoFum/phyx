#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <getopt.h>

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "bd_sim.h"
#include "log.h"

void print_help () {
    std::cout << "Birth-death simulator." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxbdsim [OPTION]... " << std::endl;
    std::cout << std::endl;
    std::cout << " -e, --extant=INT    number of extant species, alt to time" << std::endl;
    std::cout << " -t, --time=INT      depth of the tree, alt to extant" << std::endl;
    std::cout << " -b, --birth=DOUBLE  birth rate, default=1" << std::endl;
    std::cout << " -d, --death=DOUBLE  death rate, default=0" << std::endl;
    std::cout << " -n, --nreps=INT     number of replicates, default=1" << std::endl;
    std::cout << " -o, --outf=FILE     output file, stout otherwise" << std::endl;
    std::cout << " -s, --showd         show dead taxa" << std::endl;
    std::cout << " -x, --seed=INT      random number seed, clock otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

 std::string versionline("pxbdsim 0.1\nCopyright (C) 2013 FePhyFoFum\nLicense GPLv3\nwritten by Stephen A. Smith (blackrim), Joseph W. Brown");

static struct option const long_options[] =
{
    {"extant", required_argument, NULL, 'e'},
    {"time", required_argument, NULL, 't'},
    {"birth", required_argument, NULL, 'b'},
    {"death", required_argument, NULL, 'd'},
    {"nreps", required_argument, NULL, 'n'},
    {"outf", required_argument, NULL, 'o'},
    {"showd", no_argument, NULL, 's'},
    {"seed", required_argument, NULL, 'x'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'V'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool timeset = false;
    bool extantset = false;
    char * outf = NULL;
    int ext = 0;
    int nreps = 1;
    double time = 0.0;
    double birth = 1.0;
    double death = 0.0;
    bool showd = false;
    
    int seed = -1;
    while (1) {
        int oi = -1;
        int c = getopt_long(argc, argv, "e:t:b:d:n:o:x:shV", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'e':
                ext = string_to_int(optarg, "-e");
                extantset = true;
                break;
            case 't':
                time = string_to_float(optarg, "-t");
                timeset = true;
                break;
            case 'b':
                birth = string_to_float(optarg, "-b");
                if (birth <= 0) {
                    std::cerr << "Birth rate must be > 0. Exiting." << std::endl;
                    exit(0);
                }
                break;
            case 'd':
                death = string_to_float(optarg, "-d");
                if (death < 0) {
                    std::cerr << "Death rate must be >= 0. Exiting." << std::endl;
                    exit(0);
                }
                break;
            case 'n':
                nreps = string_to_int(optarg, "-n");
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                break;
            case 'x':
                seed = string_to_int(optarg, "-x");
                break;
            case 's':
                showd = true;
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
    
    if (ext == 0 && time == 0) {
        std::cerr << "You have to set -e or -t. Exiting." << std::endl;
        exit(0);
    }
    if (timeset && extantset) {
        std::cerr << "Set -e or -t, not both. Exiting." << std::endl;
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
    
    BirthDeathSimulator bd(ext, time, birth, death, seed);
    for (int i = 0; i < nreps; i++) {
        Tree * bdtr = bd.make_tree(showd);
        if (bdtr->getExtantNodeCount() > 1) {
            (*poos) << bdtr->getRoot()->getNewick(true) << ";" << std::endl;
        } else {
            (*poos) << "(" << bdtr->getRoot()->getNewick(true) << ");" << std::endl;
        }
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
