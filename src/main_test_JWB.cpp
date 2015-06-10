#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <time.h> 

using namespace std;

#include "tree_reader.h"
#include "tree.h"
#include "tree_utils.h"

bool verbose = false;

int main (int argc, char * argv[]) {
    TreeReader tr;

    if (argc != 2) {
        cout << "usage: phyx_test_JWB treefile" << endl;
        exit(0);
    }
    ifstream infile(argv[1]);
    if (!infile) {
        cerr << "Could not open treefile." << endl;
        return 1;
    }
    string line;
    vector <string> treeStrings;
    while (getline(infile, line)){
        treeStrings.push_back(line);
        if (verbose) {
            cout << "Reading tree string: "  << line << endl;
        }
    }
    infile.close();
    
    double seconds = 0;

    for (int i = 0; i < (int)treeStrings.size(); i++) {
        Tree * tree = tr.readTree(treeStrings[0]);

        cout << "Rooted status: " << is_rooted << endl;
        cout << "Tree length = "  << get_tree_length(tree) << endl;
        
        time_t start = time(NULL);
        for (int j = 0; j < 1000000; j++) {
            is_ultrametric_tips (tree);
        }
        time_t stop = time(NULL);
        seconds = difftime(stop, start);
        cout << "is_ultrametric_tips took: " << seconds << " seconds. " << endl;
        
        start = stop;
        for (int j = 0; j < 1000000; j++) {
            is_ultrametric_postorder (tree);
        }
        stop = time(NULL);
        seconds = difftime(stop, start);
        cout << "is_ultrametric_postorder took: " << seconds << " seconds. " << endl;
        
        delete tree;
    }

    return EXIT_SUCCESS;
}
