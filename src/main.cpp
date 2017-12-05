/*
 * main.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>


using namespace std;

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"

int main(int argc, char * argv[]){
    TreeReader tr;   
    string test;
       //test = "((a:2,b:2):3,(c:4,d:4):1):1;";

/*
    ifstream infile("../../big_geo/final_ml.tre.cn.rr.pr.nw.pathd8.bgstates.tre");
*/
    ifstream infile(argv[1]);
    if (!infile){
    
        cerr << "Could not open file." << endl;
        return 1;
    }
    vector<string> lines;
    string line;
    while (getline(infile, line)){
        lines.push_back(line);
    }
    infile.close();

    test = lines[0];

    Tree * tree = tr.readTree(test);
    cout << tree->getNodeCount() << endl;
/*  cout << getNewickString(tree) << endl;
    cout << tree->getRoot()->getNewick(true,"number") << ";" << endl;
    string a = "c";
    tree->pruneExternalNode(tree->getExternalNode(a));
    cout << getNewickString(tree) << endl;

    StringNodeObject sno("...a node object");
    tree->getRoot()->assocObject("test",sno);
    cout << tree->getRoot()->getNewick(true,"test") << ";" << endl;

    cout << *((StringNodeObject*) (tree->getRoot()->getObject("test"))) << endl;

    VectorNodeObject<int> vno;
    vno.push_back(1);vno.push_back(2);
    tree->getRoot()->assocObject("testvno",vno);

    cout << ((VectorNodeObject<int> *) (tree->getRoot()->getObject("testvno")))->at(0) << endl;

    a = "b";
    tree->setHeightFromRootToNodes();
    cout << tree->getExternalNode(a)->getHeight() << endl;
    cout << tree->getRoot()->getHeight() << endl;
    cout << tree->getInternalNode(0)->getHeight() << endl;
    cout << tree->getInternalNode(1)->getHeight() << endl;
    //for(int i=0;i<tree->getInternalNodeCount();i++){
        //cout << tree->getInternalNode(i).getBL() << endl;
    //}
    */
    delete tree;
    return EXIT_SUCCESS;
}
