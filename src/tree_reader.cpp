#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>

#include "node.h"
#include "tree.h"
#include "tree_reader.h"
#include "utils.h"


TreeReader::TreeReader() {}


/*
 * the tree pointer coming in should just be a new Tree()
 * we should take this out as soon as we are ready to repoint
 * the existing code to the right bits
 */

// TODO: record whether edge lengths are present, store as property
Tree * TreeReader::readTree (std::string trees) {
    Tree * tree = new Tree();
    std::string pb = trees;
    unsigned int x = 0;
    char nextChar = pb.c_str()[x];
    bool start = true;
    bool keepGoing = true;
    bool in_quote = false;
    bool hasAnnotations = false;
    bool hasInternalNodeNames = false;
    char quoteType;
    Node * currNode = NULL;
    double sumEL = 0.0;
    while (keepGoing == true) {
        //std::cout << "Working on: " << nextChar << std::endl;
        if (nextChar == '(') {
            if (start == true) {
                Node * root = new Node();
                tree->setRoot(root);
                currNode = root;
                start = false;
            } else {
                if (currNode == NULL) {
                    std::cerr << "Malformed newick string. Can read until char " << x << "." << std::endl;
                    exit(1);
                }
                Node * newNode = new Node(currNode);
                currNode->addChild(*newNode);
                currNode = newNode;
            }
        } else if (nextChar == ',') {
            currNode = currNode->getParent();
        } else if (nextChar == ')') {
            // internal named node (or more likely support annotation)
            currNode = currNode->getParent();
            x++;
            nextChar = pb.c_str()[x];
            //std::cout << "Working on: " << nextChar << std::endl;
            std::string nodeName = "";
            bool goingName = true;
            in_quote = false;
            if (nextChar == ',' || nextChar == ')' || nextChar == ':'
                || nextChar == ';' || nextChar == '[') {
                goingName = false;
            } else if (nextChar == '"' || nextChar == '\'') {
                in_quote = true;
                quoteType = nextChar;
                nodeName = nodeName + nextChar;
            }
            if (!in_quote) {
                while (goingName) {
                    nodeName = nodeName + nextChar;
                    x++;
                    nextChar = pb.c_str()[x];
                    if (nextChar == ',' || nextChar == ')' || nextChar == ':'
                        || nextChar == '[' || nextChar == ';') {
                        goingName = false;
                        break;
                    }
                }
                x--;
            } else {
                x++;
                nextChar = pb.c_str()[x];
                while (goingName) {
                    nodeName = nodeName + nextChar;
                    x++;
                    nextChar = pb.c_str()[x];
                    if (nextChar == quoteType) {
                        nodeName = nodeName + nextChar;
                        if (quoteType == '"') {
                            goingName = false;
                            break;
                        } else {
                            // check for double single quotes
                            x++;
                            nextChar = pb.c_str()[x];
                            if (nextChar != quoteType) {
                                x--;
                                nextChar = pb.c_str()[x];
                                goingName = false;
                                break;
                            }
                        }
                    }
                } 
            }// work on edge
            currNode->setName(nodeName);
            //std::cout << nodeName << std::endl;
            if (nodeName.size() > 0) {
                hasInternalNodeNames = true;
            }
            if (!in_quote) {
                //x--;
            }
        } else if (nextChar == ';') {
            keepGoing = false;
        } else if (nextChar == ':') {
            // edge length
            x++;
            nextChar = pb.c_str()[x];
            std::string edgeL = "";
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
            } // work on edge
            double edd = strtod(edgeL.c_str(), NULL);
            currNode->setBL(edd);
            sumEL += edd;
            x--;
        }
        // note/annotation
        else if (nextChar == '[') {
            hasAnnotations = true;
            x++;
            nextChar = pb.c_str()[x];
            std::string note = "";
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
            // something supposed to be here?

        }
        // external named node
        else {
            Node * newNode = new Node(currNode);
            currNode->addChild(*newNode);
            currNode = newNode;
            std::string nodeName = "";
            bool goingName = true;
            in_quote = false;
            if (nextChar == '"' || nextChar == '\'') {
                in_quote = true;
                quoteType = nextChar;
                nodeName = nodeName + nextChar;
            }
            if (!in_quote) {
                while (goingName) {
                    nodeName = nodeName + nextChar;
                    x++;
                    nextChar = pb.c_str()[x];
                    if (nextChar == ',' || nextChar == ')' || nextChar == ':'
                        || nextChar == '[') {
                        goingName = false;
                        break;
                    }
                }
                x--;
            } else {
                x++;
                nextChar = pb.c_str()[x];
                while (goingName) {
                    nodeName = nodeName + nextChar;
                    x++;
                    nextChar = pb.c_str()[x];
                    if (nextChar == quoteType) {
                        nodeName = nodeName + nextChar;
                        if (quoteType == '"') {
                            goingName = false;
                            break;
                        } else {
                            // check for double single quotes
                            x++;
                            nextChar = pb.c_str()[x];
                            if (nextChar != quoteType) {
                                x--;
                                nextChar = pb.c_str()[x];
                                goingName = false;
                                break;
                            }
                        }
                    }
                } 
            }
            //std::cout << nodeName << std::endl;
            newNode->setName(nodeName);
        }
        if (x < pb.length() - 1) { // added
            x++;
        }
        nextChar = pb.c_str()[x];
    }
    bool hasEdgeLengths = (sumEL > 0.0) ? true : false;
    tree->setEdgeLengthsPresent(hasEdgeLengths);
    //std::cout << "hasAnnotations = " << hasAnnotations << std::endl;
    tree->setNodeAnnotationsPresent(hasAnnotations);
    //std::cout << "hasInternalNodeNames = " << hasInternalNodeNames << std::endl;
    tree->setNodeNamesPresent(hasInternalNodeNames);
    tree->processRoot();
    return tree;
}


// for processing strings we know are valid tree strings
Tree * read_tree_string (std::string trees) {
    TreeReader tr;
    Tree * tree = tr.readTree(trees);
    return tree;
}


/*
 * tests the filetype by checking the first string and guessing based on
 * #nexus, ( newick
 * returns in the order above, 0, 1, 666 --- no filetype recognized
 * currently this only tests for nexus and newick
 *  if it is nexus, then the nexus reader will need 
 *  to deal with translate or not
 */
int test_tree_filetype (std::string filen) {
    std::string tline;
    std::ifstream infile(filen.c_str());
    int ret = 666; // if you get 666, there is no filetype set
    while (getline(infile, tline)) {
        if (tline.size() < 1) {
            continue;
        }
    //nexus
        if (tline[0] == '#') {
            ret = 0;
            break;
        }
    //newick
        if (tline[0] == '(') {
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
 * returns in the order above, 0, 1, 666 --- no filetype recognized
 */
int test_tree_filetype_stream (std::istream& stri, std::string& retstring) {
    if (!getline(stri, retstring)) {
        std::cout << "ERROR: end of file too soon" << std::endl;
    }
    int ret = 666; // if you get 666, there is no filetype set
    if (retstring[0] == '#') {
        ret = 0;
    } else if (retstring[0] == '(') {
        ret = 1;
    }
    return ret;
}


/**
 * this will look for the translation table if it exists
 */
bool get_nexus_translation_table (std::istream& stri, std::map<std::string, std::string> * trans,
    std::string * retstring) {
    std::string line1;
    std::string del(" \t");
    std::vector<std::string> tokens;
    bool exists = false; // is there a translation table?
    bool going = true;
    bool begintrees = false;
    bool tgoing = false;
    while (going) {
        if (!getline(stri, line1)) {
            break;
        }
        trim_spaces(line1);
        if (line1.empty()) {
            continue;
        }
        //std::cout << "Working on: " << line1 << std::endl;
        (*retstring) = line1;
        std::string uc = string_to_upper(line1);
        if (uc.find("TRANSLATE") != std::string::npos) {
            tgoing = true;
            exists = true;
            //std::cout << "Found translation table!" << std::endl;
            continue;
        } else if (begintrees == true && tgoing == false) {
            //std::cout << "No translation table present!" << std::endl;
                return false;
        }
        if (uc.find("BEGIN TREES") != std::string::npos) {
            begintrees = true;
           //std::cout << "Found Begin trees!" << std::endl;
        }
        if (tgoing == true) {
            //trimspaces and split up strings
            tokens.clear();
            tokenize(line1, tokens, del);
            size_t endcheck = line1.find(";");
            if (endcheck != std::string::npos) { // semicolon present. this is the last line.
                //std::cout << "Ending translation table!" << std::endl;
                going = false;
                (*retstring) = "";
            }
            if (tokens.size() != 1) { // not trailing lone semicolon
                for (unsigned int i=0; i < tokens.size(); i++) {
                    trim_spaces(tokens[i]);
                }
                size_t found = tokens[1].find(",");
                if (found != std::string::npos) {
                    tokens[1].erase(found, 1);
                }
                if (!going) {
                    size_t found2 = tokens[1].find(";");
                    if (found2 != std::string::npos) {
                        tokens[1].erase(found2, 1);
                    }
                }
                (*trans)[tokens[0]] = tokens[1];
                //std::cout << "tokens[0] = " << tokens[0] << "; tokens[1] = " << tokens[1] << std::endl;
            }
        }
    }
    return exists;
}


// given a line that begins with [, keep reading until a line where last char is ] (i.e., end of comment)
// NOTE: returns the end comment line (so reader will need to read in next valid line)
void process_nexus_comment (std::istream& stri, std::string& tline) {
    bool done = false;
    std::string terp = tline;
    trim_spaces(terp);
    // check single-line comment
    if (terp.back() == ']') {
        //std::cout << "single-line comment!" << std::endl;
        return;
    }
    // if not, dealing with a multi-line comment
    while (!done) {
        getline(stri, terp);
        trim_spaces(terp);
        if (!terp.empty()) {
            if (terp.back() == ']') {
                //std::cout << "found end of multi-line comment" << std::endl;
                return;
            }
        }
    }
}


bool check_nexus_comment (std::string line) {
    bool comment = false;
    trim_spaces(line);
    if (line[0] == '[') {
        comment = true;
    }
    return comment;
}


/*
 * this will read the nexus file after processing translating
 * should add some error correction code here
*/
Tree * read_next_tree_from_stream_nexus (std::istream& stri, std::string& retstring,
    bool ttexists, std::map<std::string, std::string> * trans, bool * going) {
    std::string tline;
    if (retstring.size() > 0) {
        tline = retstring;
        retstring = "";
        bool done = false;
        while (!done) {
            if (check_nexus_comment(tline)) {
                //std::cout << "yikes a comment!" << std::endl;
                process_nexus_comment(stri, tline);
                getline(stri, tline);
            } else if (tline.empty()) {
                //std::cout << "empty line. keep reading!" << std::endl;
                getline(stri, tline);
            } else {
                done = true;
            }
        }
    } else {
        bool reading = true; // continue reading if lines are empty or comments
        while (reading) {
            if (!getline(stri, tline)) {
                (*going) = false;
                return NULL;
            }
            trim_spaces(tline); // important!
            if (!tline.empty()) {
                if (check_nexus_comment(tline)) {
                    process_nexus_comment(stri, tline);
                } else {
                    reading = false;
                }
            } else {
                //std::cout << "Skipping empty line" << std::endl;
            }
        }
    }
    //std::cout << "Working on: " << tline << std::endl;
    std::string uc = string_to_upper(tline);
    if (uc.find("END;") != std::string::npos) {
        (*going) = false;
        return NULL;
    }
    //vector<string> tokens;
    //string del(" \t");
    //tokenize(tline, tokens, del);
    //string tstring(tokens[tokens.size()-1]);
    
    size_t startpos = tline.find_first_of("(");
    std::string tstring = tline.substr(startpos);
    Tree * tree;
    //std::cout << tstring << std::endl;
    TreeReader tr;
    tree = tr.readTree(tstring);
    if (ttexists) {
        for (int i=0; i < tree->getExternalNodeCount(); i++) {
            tree->getExternalNode(i)->setName((*trans)[tree->getExternalNode(i)->getName()]);
        }
    }
    return tree;
}


/*
 * this is simple as each line is a tree
 */
// adding a simple check: if line is empty, assume we're done
Tree * read_next_tree_from_stream_newick (std::istream& stri, std::string& retstring, bool * going) {
    std::string tline;
    if (retstring.size() > 0) {
        tline = retstring;
        retstring = "";
    } else if (!getline(stri, tline)) {
        (*going) = false;
        return NULL;
    }
    if (tline.size() == 0) {
        //std::cout << "You've got yerself an empty line, there." << std::endl;
        (*going) = false;
        return NULL;
    }
    Tree * tree;
    TreeReader tr;
    tree = tr.readTree(tline);
    return tree;
}
