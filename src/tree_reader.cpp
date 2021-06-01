#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>

#include "node.h"
#include "tree.h"
#include "tree_reader.h"
#include "utils.h"


TreeReader::TreeReader() = default;


/*
 * the tree pointer coming in should just be a new Tree()
 * we should take this out as soon as we are ready to repoint
 * the existing code to the right bits
*/
Tree * TreeReader::readTree (const std::string& pb) {
    auto * tree = new Tree();
    //std::string pb = trees;
    unsigned int x = 0;
    char nextChar = pb.c_str()[x];
    bool start = true;
    bool keepGoing = true;
    bool in_quote = false;
    bool hasAnnotations = false;
    bool hasInternalNodeNames = false;
    char quoteType = '\0';
    Node * currNode = nullptr;
    double sumEL = 0.0;
    while (keepGoing) {
        //std::cout << "Working on: " << nextChar << std::endl;
        if (nextChar == '(') {
            if (start) {
                auto * root = new Node();
                tree->setRoot(root);
                currNode = root;
                start = false;
            } else {
                if (currNode == nullptr) {
                    std::cerr << "Malformed newick string. Can read until char " << x << "." << std::endl;
                    exit(1);
                }
                auto * newNode = new Node(currNode);
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
            std::string nodeName;
            bool goingName = true;
            in_quote = false;
            if (nextChar == ',' || nextChar == ')' || nextChar == ':'
                || nextChar == ';' || nextChar == '[') {
                goingName = false;
            } else if (nextChar == '"' || nextChar == '\'') {
                in_quote = true;
                quoteType = nextChar;
                nodeName += nextChar;
            }
            if (!in_quote) {
                while (goingName) {
                    nodeName += nextChar;
                    x++;
                    nextChar = pb.c_str()[x];
                    if (nextChar == ',' || nextChar == ')' || nextChar == ':'
                        || nextChar == '[' || nextChar == ';') {
                        break;
                    }
                }
                x--;
            } else {
                x++;
                nextChar = pb.c_str()[x];
                while (goingName) {
                    nodeName += nextChar;
                    x++;
                    nextChar = pb.c_str()[x];
                    if (nextChar == quoteType) {
                        nodeName += nextChar;
                        if (quoteType == '"') {
                            break;
                        }
                        // check for double single quotes
                        x++;
                        nextChar = pb.c_str()[x];
                        if (nextChar != quoteType) {
                            x--;
                            break;
                        }
                    }
                } 
            }// work on edge
            currNode->setName(nodeName);
            //std::cout << nodeName << std::endl;
            if (!nodeName.empty()) {
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
            std::string edgeL;
            bool goingName = true;
            while (goingName) {
                edgeL += nextChar;
                x++;
                nextChar = pb.c_str()[x];
                if (nextChar == ',' || nextChar == ')' || nextChar == ':'
                    || nextChar == ';'|| nextChar == '[') {
                    break;
                }
            } // work on edge
            double edd = strtod(edgeL.c_str(), nullptr);
            currNode->setBL(edd);
            sumEL += edd;
            x--;
        }
        // note/annotation
        else if (nextChar == '[') {
            hasAnnotations = true;
            x++;
            nextChar = pb.c_str()[x];
            std::string note;
            bool goingNote = true;
            while (goingNote) {
                note += nextChar;
                x++;
                nextChar = pb.c_str()[x];
                if (nextChar == ']' ) {
                    break;
                }
            }
            currNode->setComment(note);
        } else if (nextChar == ' ') {
            // just skip internal whitespace
        }
        // external named node
        else {
            auto * newNode = new Node(currNode);
            currNode->addChild(*newNode);
            currNode = newNode;
            std::string nodeName;
            bool goingName = true;
            in_quote = false;
            if (nextChar == '"' || nextChar == '\'') {
                in_quote = true;
                quoteType = nextChar;
                nodeName += nextChar;
            }
            if (!in_quote) {
                while (goingName) {
                    nodeName += nextChar;
                    x++;
                    nextChar = pb.c_str()[x];
                    if (nextChar == ',' || nextChar == ')' || nextChar == ':'
                        || nextChar == '[') {
                        break;
                    }
                }
                x--;
            } else {
                x++;
                nextChar = pb.c_str()[x];
                while (goingName) {
                    nodeName += nextChar;
                    x++;
                    nextChar = pb.c_str()[x];
                    if (nextChar == quoteType) {
                        nodeName += nextChar;
                        if (quoteType == '"') {
                            break;
                        }
                        // check for double single quotes
                        x++;
                        nextChar = pb.c_str()[x];
                        if (nextChar != quoteType) {
                            x--;
                            break;
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
    //bool hasEdgeLengths = (sumEL > 0.0) ? true : false;
    bool hasEdgeLengths = sumEL > 0.0;
    tree->setEdgeLengthsPresent(hasEdgeLengths);
    //std::cout << "hasAnnotations = " << hasAnnotations << std::endl;
    tree->setNodeAnnotationsPresent(hasAnnotations);
    //std::cout << "hasInternalNodeNames = " << hasInternalNodeNames << std::endl;
    tree->setNodeNamesPresent(hasInternalNodeNames);
    tree->processRoot();
    return tree;
}


// for processing strings we know are valid tree strings
Tree * read_tree_string (std::string pb) {
    TreeReader tr;
    Tree * tree = tr.readTree(std::move(pb));
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
// hrm this does not seem to be used anymore
int test_tree_filetype (const std::string& filen) {
    std::string tline;
    std::ifstream infile(filen.c_str());
    int ret = 666; // if you get 666, there is no filetype set
    while (getline_safe(infile, tline)) {
        if (tline.empty()) {
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
    if (!getline_safe(stri, retstring)) {
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
        if (!getline_safe(stri, line1)) {
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
        }
        if (begintrees && !tgoing) {
            //std::cout << "No translation table present!" << std::endl;
            return false;
        }
        if (uc.find("BEGIN TREES") != std::string::npos) {
            begintrees = true;
           //std::cout << "Found Begin trees!" << std::endl;
        }
        if (tgoing) {
            //trimspaces and split up strings
            tokens.clear();
            tokenize(line1, tokens, del);
            size_t endcheck = line1.find(';');
            if (endcheck != std::string::npos) { // semicolon present. this is the last line.
                //std::cout << "Ending translation table!" << std::endl;
                going = false;
                (*retstring) = "";
            }
            if (tokens.size() != 1) { // not trailing lone semicolon
                for (auto & tk : tokens) {
                    trim_spaces(tk);
                }
                size_t found = tokens[1].find(',');
                if (found != std::string::npos) {
                    tokens[1].erase(found, 1);
                }
                if (!going) {
                    size_t found2 = tokens[1].find(';');
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


/*
 * this will read the nexus file after processing translating
 * should add some error correction code here
*/
Tree * read_next_tree_from_stream_nexus (std::istream& stri, std::string& retstring,
    bool ttexists, std::map<std::string, std::string> * trans, bool * going) {
    std::string tline;
    if (!retstring.empty()) {
        tline = retstring;
        retstring = "";
        bool done = false;
        while (!done) {
            if (check_nexus_comment(tline)) {
                //std::cout << "yikes a comment!" << std::endl;
                process_nexus_comment(stri, tline);
                getline_safe(stri, tline);
            } else if (tline.empty()) {
                //std::cout << "empty line. keep reading!" << std::endl;
                getline_safe(stri, tline);
            } else {
                done = true;
            }
        }
    } else {
        bool reading = true; // continue reading if lines are empty or comments
        while (reading) {
            if (!getline_safe(stri, tline)) {
                (*going) = false;
                return nullptr;
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
        return nullptr;
    }
    //vector<string> tokens;
    //string del(" \t");
    //tokenize(tline, tokens, del);
    //string tstring(tokens[tokens.size()-1]);
    
    size_t startpos = tline.find_first_of('(');
    std::string tstring = tline.substr(startpos);
    Tree * tree;
    //std::cout << tstring << std::endl;
    TreeReader tr;
    tree = tr.readTree(tstring);
    if (ttexists) {
        for (int i = 0; i < tree->getExternalNodeCount(); i++) {
            tree->getExternalNode(i)->setName((*trans)[tree->getExternalNode(i)->getName()]);
        }
    }
    return tree;
}


/*
 * this is simple as each line is a tree
 */
// adding a simple check: if line is empty, assume we're done
// TODO: deal with trees with internal line breaks (phylip does this?)
Tree * read_next_tree_from_stream_newick (std::istream& stri, std::string& retstring, bool * going) {
    std::string tline;
    if (!retstring.empty()) {
        tline = retstring;
        retstring = "";
    } else if (!getline_safe(stri, tline)) {
        (*going) = false;
        return nullptr;
    }
    trim_spaces(tline);
    
    // hrm do we want to allow empty lines in between trees? i think so?
    if (tline.empty()) {
        //std::cout << "You've got yerself an empty line, there." << std::endl;
        bool done = false;
        while (!done) {
            if (!getline_safe(stri, tline)) {
                (*going) = false;
                return nullptr;
            }
            trim_spaces(tline);
            if (!tline.empty()) {
                done = true;
            }
        }
    }
    
    if (tline.back() != ';') {
        bool done = false;
        std::string terp;
        while (!done) {
            if (!getline_safe(stri, terp)) {
                std::cerr << "Error: malformed tree string (missing trailing semicolon). Exiting." << std::endl;
                exit(1);
            } else {
                trim_spaces(terp);
                tline += terp;
                if (tline.back() == ';') {
                    done = true;
                }
            }
        }
    }
    
    Tree * tree;
    TreeReader tr;
    tree = tr.readTree(tline);
    return tree;
}
