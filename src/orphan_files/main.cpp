#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"

int main (int argc, char * argv[]) {
    TreeReader tr;   
    std::string test;
       //test = "((a:2,b:2):3,(c:4,d:4):1):1;";

/*
    std::ifstream infile("../../big_geo/final_ml.tre.cn.rr.pr.nw.pathd8.bgstates.tre");
*/
    std::ifstream infile(argv[1]);
    if (!infile){
        std::cerr << "Could not open file." << std::endl;
        return 1;
    }
    std::vector<std::string> lines;
    std::string line;
    while (getline(infile, line)){
        lines.push_back(line);
    }
    infile.close();

    test = lines[0];

    Tree * tree = tr.readTree(test);
    std::cout << tree->getNodeCount() << std::endl;
/*  std::cout << getNewickString(tree) << std::endl;
    std::cout << tree->getRoot()->getNewick(true,"number") << ";" << std::endl;
    std::string a = "c";
    tree->pruneExternalNode(tree->getExternalNode(a));
    std::cout << getNewickString(tree) << std::endl;

    StringNodeObject sno("...a node object");
    tree->getRoot()->assocObject("test",sno);
    std::cout << tree->getRoot()->getNewick(true,"test") << ";" << std::endl;

    std::cout << *((StringNodeObject*) (tree->getRoot()->getObject("test"))) << std::endl;

    VectorNodeObject<int> vno;
    vno.push_back(1);vno.push_back(2);
    tree->getRoot()->assocObject("testvno",vno);

    std::cout << ((VectorNodeObject<int> *) (tree->getRoot()->getObject("testvno")))->at(0) << std::endl;

    a = "b";
    tree->setHeightFromRootToNodes();
    std::cout << tree->getExternalNode(a)->getHeight() << std::endl;
    std::cout << tree->getRoot()->getHeight() << std::endl;
    std::cout << tree->getInternalNode(0)->getHeight() << std::endl;
    std::cout << tree->getInternalNode(1)->getHeight() << std::endl;
    //for(int i=0;i<tree->getInternalNodeCount();i++){
        //cout << tree->getInternalNode(i).getBL() << std::endl;
    //}
    */
    delete tree;
    return EXIT_SUCCESS;
}
