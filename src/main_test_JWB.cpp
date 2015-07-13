#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <time.h> 

using namespace std;

#include "utils.h"
#include "tree_reader.h"
#include "tree.h"
#include "tree_utils.h"

bool verbose = false;

//g++ -std=c++11 main_test_JWB.cpp utils.cpp tree_utils.cpp tree.cpp tree_reader.cpp superdouble.cpp node.cpp -o test_JWB

int main (int argc, char * argv[]) {
    
    // testing labels
    string test = "Anomalurus_sp__GP_2005_scaly-tailed_squirrel";
    string fixed = get_valid_nexus_label(test);
    string test2 = "Anomalurus_sp__GP_2005_John's_scaly-tailed_squirrel";
    string fixed2 = get_valid_nexus_label(test2);
    
    string test3 = "Anomalurus_sp__GP_2005_scaly-tailed_squirrel+";
    string fixed3 = get_valid_nexus_label(test3);
    string fixed4 = get_valid_newick_label(test3);
    
    // no invalid characters, but spaces
    string test5 = "Anomalurus sp  GP 2005 scaly tailed squirrel";
    string fixed5 = get_valid_nexus_label(test5);
    
    cout << "Original:  " << test << endl;
    cout << "Fixed:     " << fixed << endl;
    cout << "Original2: " << test2 << endl;
    cout << "Fixed2:    " << fixed2 << endl;
    
    cout << "Original3: " << test3 << endl;
    cout << "Fixed3:    " << fixed3 << endl;
    cout << "Fixed4:    " << fixed4 << endl;
    
    cout << "Original5: " << test5 << endl;
    cout << "Fixed5:    " << fixed5 << endl;
    
    if (false) {
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

        for (unsigned int i = 0; i < treeStrings.size(); i++) {
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
    }

    return EXIT_SUCCESS;
}
