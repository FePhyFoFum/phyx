

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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

    for (int i = 0; i < (int)treeStrings.size(); i++) {
        Tree * tree = tr.readTree(treeStrings[0]);

        cout << "Rooted status: " << is_rooted << endl;
        cout << "Tree length = "  << get_tree_length(tree) << endl;
        if (is_ultrametric_tips (tree)) {
            cout << "Tree is ultrametric (tip-to-root paths)!" << endl;
        } else {
            cout << "Tree is NOT ultrametric (tip-to-root paths). Wah-wha..." << endl;
        }
        if (is_ultrametric_postorder (tree)) {
            cout << "Tree is ultrametric (postorder)!" << endl;
        } else {
            cout << "Tree is NOT ultrametric (postorder). Wah-wha..." << endl;
        }
        delete tree;
    }

	return EXIT_SUCCESS;
}
