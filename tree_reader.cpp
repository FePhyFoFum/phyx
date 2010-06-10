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

using namespace std;

#include "node.h"
#include "tree.h"
#include "tree_reader.h"

TreeReader::TreeReader(){}

/*
 * the tree pointer coming in should just be a new Tree()
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
