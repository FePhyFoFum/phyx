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
#include "citations.h"


void print_help ();
std::string get_version_line ();

void print_help () {
    std::cout << std::endl;
    std::cout << "Birth-death tree simulator." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: pxbdsim [OPTIONS]..." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -e, --extant=INT    number of extant species, alt to time" << std::endl;
    std::cout << " -t, --time=DOUBLE   timespan of simulation (age of root), alt to extant" << std::endl;
    std::cout << " -b, --birth=DOUBLE  birth rate, default=1" << std::endl;
    std::cout << " -d, --death=DOUBLE  death rate, default=0" << std::endl;
    std::cout << " -n, --nreps=INT     number of replicates, default=1" << std::endl;
    std::cout << " -s, --showextinct   show lineages that went extinct, default=false" << std::endl;
    std::cout << " -x, --seed=INT      random number seed, clock otherwise" << std::endl;
    std::cout << " -v, --verbose       print per-tree simulation summary (to cerr)" << std::endl;
    std::cout << " -o, --outf=FILE     output file, STOUT otherwise" << std::endl;
    std::cout << " -h, --help          display this help and exit" << std::endl;
    std::cout << " -V, --version       display version and exit" << std::endl;
    std::cout << " -C, --citation      display phyx citation and exit" << std::endl;
    std::cout << std::endl;
    std::cout << "Report bugs to: <https://github.com/FePhyFoFum/phyx/issues>" << std::endl;
    std::cout << "phyx home page: <https://github.com/FePhyFoFum/phyx>" << std::endl;
}

 std::string get_version_line () {
    std::string vl = "pxbdsim 1.3.2\n";
    vl += "Copyright (C) 2013-2025 FePhyFoFum\n";
    vl += "License GPLv3\n";
    vl += "Written by Stephen A. Smith (blackrim), Joseph W. Brown";
    return vl;
}

static struct option const long_options[] =
{
    {"extant", required_argument, nullptr, 'e'},
    {"time", required_argument, nullptr, 't'},
    {"birth", required_argument, nullptr, 'b'},
    {"death", required_argument, nullptr, 'd'},
    {"nreps", required_argument, nullptr, 'n'},
    {"outf", required_argument, nullptr, 'o'},
    {"showextinct", no_argument, nullptr, 's'},
    {"seed", required_argument, nullptr, 'x'},
    {"verbose", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {"version", no_argument, nullptr, 'V'},
    {"citation", no_argument, nullptr, 'C'},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char * argv[]) {
    
    log_call(argc, argv);
    
    bool outfileset = false;
    bool timeset = false;
    bool extantset = false;
    char * outf = nullptr;
    int ext = 0;
    int nreps = 1;
    double time = 0.0;
    double birth = 1.0;
    double death = 0.0;
    bool showd = false;
    bool verbose = false;
    long int seed = -1;
    bool argspresent = false;
    
    while (true) {
        int oi = -1;
        int c = getopt_long(argc, argv, "e:t:b:d:n:o:x:vshVC", long_options, &oi);
        if (c == -1) {
            break;
        }
        switch(c) {
            case 'e':
                ext = string_to_int(optarg, "-e");
                extantset = true;
                argspresent = true;
                break;
            case 't':
                time = string_to_double(optarg, "-t");
                timeset = true;
                argspresent = true;
                break;
            case 'b':
                birth = string_to_double(optarg, "-b");
                if (birth <= 0) {
                    std::cerr << "Error: birth rate must be > 0. Exiting." << std::endl;
                    exit(0);
                }
                argspresent = true;
                break;
            case 'd':
                death = string_to_double(optarg, "-d");
                if (death < 0) {
                    std::cerr << "Error: death rate must be >= 0. Exiting." << std::endl;
                    exit(0);
                }
                argspresent = true;
                break;
            case 'n':
                nreps = string_to_int(optarg, "-n");
                argspresent = true;
                break;
            case 'v':
                verbose = true;
                argspresent = true;
                break;
            case 'o':
                outfileset = true;
                outf = strdup(optarg);
                argspresent = true;
                break;
            case 'x':
                seed = string_to_long_int(optarg, "-x");
                argspresent = true;
                break;
            case 's':
                showd = true;
                argspresent = true;
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
    
    if (!argspresent) {
        print_help();
        exit(1);
    }
    
    if (!extantset && !timeset) {
        std::cerr << "Error: you have to set -e or -t. Exiting." << std::endl;
        exit(0);
    }
    if (timeset && extantset) {
        std::cerr << "Error: set -e or -t, not both. Exiting." << std::endl;
        exit(0);
    }
    
     std::ostream * poos = nullptr;
     std::ofstream * ofstr = nullptr;
    
    if (outfileset) {
        ofstr = new std::ofstream(outf);
        poos = ofstr;
    } else {
        poos = &std::cout;
    }
    
    BirthDeathSimulator bd(ext, time, birth, death, seed);
    for (int i = 0; i < nreps; i++) {
        Tree * bdtr = bd.make_tree(showd);
        
        // only extinct-pruned scenarios can produce a single terminal tree
        // tree writer has a problem with this, so need to manually fix
        if (!showd && bdtr->getExtantNodeCount() == 1) {
            (*poos) << "(" << bdtr->getRoot()->getNewick(true) << ");" << std::endl;
        } else {
            (*poos) << bdtr->getRoot()->getNewick(true) << ";" << std::endl;
        }
        // print simulation summary
        if (verbose) {
            std::cerr << bd.get_sim_summary() << std::endl;
        }
    }
    if (outfileset) {
        ofstr->close();
        delete poos;
    }
    return EXIT_SUCCESS;
}
