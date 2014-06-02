/*
 * tree_reader.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>

using namespace std;

#include "node.h"
#include "tree.h"
#include "tree_reader.h"
#include "utils.h"

TreeReader::TreeReader(){}

/*
 * the tree pointer coming in should just be a new Tree()
 * we should take this out as soon as we are ready to repoint
 * the existing code to teh right bits
 */
Tree * TreeReader::readTree(string trees){
    Tree * tree = new Tree();
    string pb = trees;
    unsigned int x = 0;
    char nextChar = pb.c_str()[x];
    bool start = true;
    bool keepGoing = true;
    Node * currNode = NULL;
    while (keepGoing == true) {
	if (nextChar == '(') {
	    if (start == true) {
		Node * root = new Node();
		tree->setRoot(root);
		currNode = root;
		start = false;
	    } else {
		Node * newNode = new Node(currNode);
		currNode->addChild(*newNode);
		currNode = newNode;
	    }
	} else if (nextChar == ',') {
	    currNode = currNode->getParent();
	} else if (nextChar == ')') {
	    currNode = currNode->getParent();
	    x++;
	    nextChar = pb.c_str()[x];
	    string nam = "";
	    bool goingName = true;
	    if (nextChar == ',' || nextChar == ')' || nextChar == ':'
		|| nextChar == ';'|| nextChar == '[') {
		goingName = false;
	    }
	    while (goingName == true) {
		nam = nam + nextChar;
		x++;
		nextChar = pb.c_str()[x];
		if (nextChar == ',' || nextChar == ')' || nextChar == ':'
		    || nextChar == ';'|| nextChar == '[') {
		    goingName = false;
		    break;
		}
	    }// work on edge
	    currNode->setName(nam);
	    x--;
	} else if (nextChar == ';') {
	    keepGoing = false;
	} else if (nextChar == ':') {
	    x++;
	    nextChar = pb.c_str()[x];
	    string edgeL = "";
	    bool goingName = true;
	    while (goingName == true) {
		edgeL = edgeL + nextChar;
		x++;
		nextChar = pb.c_str()[x];
		if (nextChar == ',' || nextChar == ')' || nextChar == ':'
		    || nextChar == ';'|| nextChar == '[') {
		    goingName = false;
		    break;
		}
	    }// work on edge
	    double edd = strtod(edgeL.c_str(),NULL);
	    currNode->setBL(edd);
	    x--;
	}
	//note
	else if (nextChar == '[') {
	    x++;
	    nextChar = pb.c_str()[x];
	    string note = "";
	    bool goingNote = true;
	    while (goingNote == true) {
		note = note + nextChar;
		x++;
		nextChar = pb.c_str()[x];
		if (nextChar == ']' ) {
		    goingNote = false;
		    break;
		}
	    }
	    currNode->setComment(note);
	} else if (nextChar == ' ') {
	    
	}
	// external named node
	else {
	    Node * newNode = new Node(currNode);
	    currNode->addChild(*newNode);
	    currNode = newNode;
	    string nodeName = "";
	    bool goingName = true;
	    while (goingName == true) {
		nodeName = nodeName + nextChar;
		x++;
		nextChar = pb.c_str()[x];
		if (nextChar == ',' || nextChar == ')' || nextChar == ':' || nextChar == '[') {
		    goingName = false;
		    break;
		}
	    }
	    newNode->setName(nodeName);
	    x--;
	}
	if (x < pb.length() - 1)//added
	    x++;
	nextChar = pb.c_str()[x];
    }
    tree->processRoot();
    return tree;
}

Tree * read_tree_string(string trees){
    Tree * tree = new Tree();
    string pb = trees;
    unsigned int x = 0;
    char nextChar = pb.c_str()[x];
    bool start = true;
    bool keepGoing = true;
    Node * currNode = NULL;
    while (keepGoing == true) {
	if (nextChar == '(') {
	    if (start == true) {
		Node * root = new Node();
		tree->setRoot(root);
		currNode = root;
		start = false;
	    } else {
		Node * newNode = new Node(currNode);
		currNode->addChild(*newNode);
		currNode = newNode;
	    }
	} else if (nextChar == ',') {
	    currNode = currNode->getParent();
	} else if (nextChar == ')') {
	    currNode = currNode->getParent();
	    x++;
	    nextChar = pb.c_str()[x];
	    string nam = "";
	    bool goingName = true;
	    if (nextChar == ',' || nextChar == ')' || nextChar == ':'
		|| nextChar == ';'|| nextChar == '[') {
		goingName = false;
	    }
	    while (goingName == true) {
		nam = nam + nextChar;
		x++;
		nextChar = pb.c_str()[x];
		if (nextChar == ',' || nextChar == ')' || nextChar == ':'
		    || nextChar == ';'|| nextChar == '[') {
		    goingName = false;
		    break;
		}
	    }// work on edge
	    currNode->setName(nam);
	    x--;
	} else if (nextChar == ';') {
	    keepGoing = false;
	} else if (nextChar == ':') {
	    x++;
	    nextChar = pb.c_str()[x];
	    string edgeL = "";
	    bool goingName = true;
	    while (goingName == true) {
		edgeL = edgeL + nextChar;
		x++;
		nextChar = pb.c_str()[x];
		if (nextChar == ',' || nextChar == ')' || nextChar == ':'
		    || nextChar == ';'|| nextChar == '[') {
		    goingName = false;
		    break;
		}
	    }// work on edge
	    double edd = strtod(edgeL.c_str(),NULL);
	    currNode->setBL(edd);
	    x--;
	}
	//note
	else if (nextChar == '[') {
	    x++;
	    nextChar = pb.c_str()[x];
	    string note = "";
	    bool goingNote = true;
	    while (goingNote == true) {
		note = note + nextChar;
		x++;
		nextChar = pb.c_str()[x];
		if (nextChar == ']' ) {
		    goingNote = false;
		    break;
		}
	    }
	    currNode->setComment(note);
	} else if (nextChar == ' ') {
	    
	}
	// external named node
	else {
	    Node * newNode = new Node(currNode);
	    currNode->addChild(*newNode);
	    currNode = newNode;
	    string nodeName = "";
	    bool goingName = true;
	    while (goingName == true) {
		nodeName = nodeName + nextChar;
		x++;
		nextChar = pb.c_str()[x];
		if (nextChar == ',' || nextChar == ')' || nextChar == ':' || nextChar == '[') {
		    goingName = false;
		    break;
		}
	    }
	    newNode->setName(nodeName);
	    x--;
	}
	if (x < pb.length() - 1)//added
	    x++;
	nextChar = pb.c_str()[x];
    }
    tree->processRoot();
    return tree;
}

/*
 * tests the filetype by checking the first string and guessing based on
 * #nexus, ( newick
 * returns in the order above, 0 ,1 , 666 --- no filetype recognized
 * currently this only tests for nexus and newick
 *  if it is nexus, then the nexus reader will need 
 *  to deal with translate or not
 */
int test_tree_filetype(string filen){
    string tline;
    ifstream infile(filen.c_str());
    int ret = 666; // if you get 666, there is no filetype set
    while (getline(infile,tline)){
        if (tline.size() < 1){
            continue;
        }
	//nexus
        if (tline[0] == '#'){
            ret = 0;
            break;
        }
	//newick
        if (tline[0] == '('){
            ret = 1;
            break;
        }
        break;
    }
    infile.close();
    return ret;
}

/* tests the filetype by checking the first string and guessing based on 
 * #nexus, ( newick
 * returns in the order above, 0 ,1, 666 --- no filetype recognized
 */
int test_tree_filetype_stream(istream & stri, string & retstring){
    if (!getline(stri, retstring)){
        cout << "ERROR: end of file too soon" << endl;
    }
    int ret = 666; // if you get 666, there is no filetype set
    if (retstring[0] == '#'){
        ret = 0;
    }else if (retstring[0] == '('){
        ret = 1;
    }
    return ret;
}

/**
 * this will look for the translation table i it exists
 * 
 */
bool get_nexus_translation_table(istream & stri, map<string, string> * trans, vector<string> * retstrings){
    string line1;
    bool exists = false;
    bool going = true;
    bool begintrees = false;
    bool tgoing = false;
    while (going){
	if(!getline(stri,line1)){
	    break;
	}
	(*retstrings).push_back(line1);
	string uc(line1);
	transform(uc.begin(), uc.end(), uc.begin(), ::toupper);
	if(uc.find("TRANSLATE") != string::npos){
	    tgoing = true;
	    exists = true;
	    continue;
	}else if(begintrees == true && tgoing == false){
	    return false;
	}
	if(uc.find("BEGIN TREES") != string::npos){
	    begintrees = true;
	}
	if(tgoing == true){
	    //trimspaces and split up strings
	    vector<string> tokens;
	    string del(" \t");
	    tokens.clear();
	    tokenize(line1,tokens,del);
	    for(int i=0;i<tokens.size();i++){
		trim_spaces(tokens[i]);
	    }
	    size_t found = tokens[1].find(",");
	    if (found != string::npos)
		tokens[1].erase(found,1);
	    (*trans)[tokens[0]] = tokens[1];
	    size_t found2 = tokens[1].find(";");
	    if(found2 != string::npos){
		tokens[1].erase(found2,1);
		(*trans)[tokens[0]] = tokens[1];
		going = false;
	    }
	}
    }
    return exists;
}


/*
 * this will read the nexus file after processing translating
 * should add some error correction code here
 */

Tree * read_next_tree_from_stream_nexus(istream & stri, string & retstring, bool ttexists, map<string,string> * trans, bool * going){
     string tline;
     if(!getline(stri,tline)){
	 (*going) = false;
	 return NULL;
     }
     string uc(tline);
     transform(uc.begin(),uc.end(),uc.begin(),::toupper);
     if(uc.find("END;")!=string::npos){
	 (*going) = false;
	 return NULL;
     }
     vector<string> tokens;
     string del(" \t");
     tokenize(tline,tokens,del);
     string tstring(tokens[tokens.size()-1]);
     Tree * tree;
     //cout << tstring << endl;
     tree = read_tree_string(tstring);
     if(ttexists){
	 for(int i=0;i<tree->getExternalNodeCount();i++){
	     tree->getExternalNode(i)->setName((*trans)[tree->getExternalNode(i)->getName()]);
	 }
     }
     return tree;
}

/*
 * this is simple as each line is a tree
 *
 */
Tree * read_next_tree_from_stream_newick(istream & stri, string & retstring, bool * going){
    string tline;
    if(retstring.size() > 0){
	tline = retstring;
	retstring = "";
    }else if(!getline(stri,tline)){
	(*going) = false;
	return NULL;
    }
    Tree * tree;
    tree = read_tree_string(tline);
    return tree;
}
