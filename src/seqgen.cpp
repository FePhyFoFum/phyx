/*
 * seqgen.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: joe
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <map>
#include <iterator>
#include <cstring>
#include <getopt.h>
#include <stdio.h>
#include <nlopt.hpp>
#include <math.h>
#include <armadillo>
#include <random>


using namespace arma;
using namespace std;

#include "seqgen.h"
#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"
#include "node.h"
#include "tree.h"
#include "tree_reader.h"

/*Default Rate Matrix looks like this, I don't know why but I always go A,T,C,G
 *
 *    A   T   C   G
 * A -1  .33 .33 .33
 * T .33 -1  .33 .33
 * G .33 .33  -1 .33
 * G .33 .33 .33  -1
 */


/*Use the P matrix probabilities and randomly draw numbers to see
 *if each individual state will undergo some type of change
 *
 *
 */
void SeqSim(string& Ancestor, vector< vector<double> > Matrix){


    string hold = ""; //for the stupid to double thing I always do
    string newstring = "";
    for (int i = 0; i < Ancestor.size(); i++){
        float RandNumb = 0.0;
        float ChanceA = 0.0;
        float ChanceT = 0.0;
        float ChanceC = 0.0;
        float ChanceG = 0.0;
        random_device rand_dev;
        mt19937 generator(rand_dev());
        std::uniform_int_distribution<int>  distr(0, 100000);
    	hold = to_string(distr(generator));
    	RandNumb = stod(hold);
    	RandNumb /= 100000;
    	if (Ancestor[i] == 'A'){
    		ChanceA = abs(Matrix[0][0]);
    		ChanceT = Matrix[0][1] /abs(Matrix[0][0]);
    		ChanceC = (Matrix[0][2] /abs(Matrix[0][0]) + ChanceT);
    		ChanceG = (Matrix[0][3] /abs(Matrix[0][0]) + ChanceC);
    		if(RandNumb < ChanceT){
    			newstring += "T";
    		}else if(RandNumb > ChanceT && RandNumb < ChanceC){

    			newstring += "C";
    		}else if(RandNumb > ChanceC && RandNumb < ChanceG){

    			newstring += "G";
    		}else{

    			newstring += "A";
    		}
    	}else if (Ancestor[i] == 'T'){

    		ChanceT = abs(Matrix[1][1]);
    		ChanceA = Matrix[1][0] / ChanceT;
    		ChanceC = (Matrix[1][2] / ChanceT + ChanceA);
    		ChanceG = (Matrix[1][3] / ChanceT + ChanceC);
    		if(RandNumb < ChanceA){
    			newstring += "A";
    		}else if(RandNumb > ChanceA && RandNumb < ChanceC){
    			newstring += "C";
    		}else if(RandNumb > ChanceC && RandNumb < ChanceG){

    			newstring += "G";
    		}else{

    			newstring += "T";
    		}
    	}else if (Ancestor[i] == 'C'){

    		ChanceC = abs(Matrix[2][2]);
    		ChanceA = Matrix[2][0] / ChanceC;
    		ChanceT = (Matrix[2][1] /ChanceC + ChanceA);
    		ChanceG = (Matrix[2][3] /ChanceC + ChanceT);
    		if(RandNumb < ChanceA){
    			newstring += "A";
    		}else if(RandNumb > ChanceA && RandNumb < ChanceT){

    			newstring += "T";
    		}else if(RandNumb > ChanceT && RandNumb < ChanceG){

    			newstring += "G";
    		}else{

    			newstring += "C";
    		}
    	}else{

    		ChanceG = abs(Matrix[3][3]);
    		ChanceA = Matrix[3][0] / ChanceG;
    		ChanceT = (Matrix[3][1] / ChanceG + ChanceA);
    		ChanceC = (Matrix[3][2] / ChanceG + ChanceC);
    		if(RandNumb < ChanceA){
    			newstring += "A";
    		}else if(RandNumb > ChanceA && RandNumb < ChanceT){

    			newstring += "T";
    		}else if(RandNumb > ChanceT && RandNumb < ChanceC){

    			newstring += "C";
    		}else{

    			newstring += "G";
    		}
    	}

    	//cout << (RandNumb / 10000) << endl;
    }
    Ancestor = newstring;

}
/*
 * Calculate the Q Matrix (Substitution rate matrix)
 */
vector< vector<double> > CalQ(vector< vector<double> >& rmatrix, vector<double>& basefreq){

	vector< vector<double> > bigpi(4, vector<double>(4, 1.0));
	vector< vector<double> > t(4, vector<double>(4, 0.0));
	double tscale = 0.0;
    for (int i = 0; i < rmatrix.size(); i++){
    	for (int j = 0; j < rmatrix.size(); j++){

    		if (i != j){

    			bigpi[i][j] *= basefreq[i] * basefreq[j] * rmatrix[i][j];
    			tscale += bigpi[i][j];

    		}else{
    			bigpi[i][j] = 0.0;
    		}
    	}
    }
    for (int i = 0; i < rmatrix.size(); i++){
    	for (int j = 0; j < rmatrix.size(); j++){

    		if (i != j){

    			bigpi[i][j] /= tscale;

    		}else{
    			//set the diagnols to zero
    			bigpi[i][j] = 0.0;
    		}
    	}
    }
    for (int i = 0; i < rmatrix.size(); i++){
    	double diag = 0.0;
    	for (int j = 0; j < rmatrix.size(); j++){

    		if (i != j){

    			diag -= bigpi[i][j];

    		}
    	}
    	bigpi[i][i] = diag;
    }
    //Divide and Transpose
    for (int i = 0; i < rmatrix.size(); i++){
    	for (int j = 0; j < rmatrix.size(); j++){

    		bigpi[i][j] /= basefreq[i];
    	}
    }

	return bigpi;
}
/*Calculate the P Matrix (Probability Matrix
 * Changes to armadillos format then back I don't like the way could be more
 * efficient but yeah...
 */
vector< vector<double> > PCalq(vector< vector<double> > QMatrix, float br){

	vector< vector<double> > Pmatrix(4, vector<double>(4, 0.0));
	mat A = randn<mat>(4,4);
	mat B = randn<mat>(4,4);
	int count = 0;
	//Q * t moved into Matrix form for armadillo
   for (int i = 0; i < QMatrix.size(); i++){
    	for (int j = 0; j < QMatrix.size(); j++){

    		A[count] = (QMatrix[i][j] * br);
    		count++;
    	}
    }
   //exponentiate the matrix
   B = expmat(A);
   //cout << B << endl;
   count = 0;
   //convert the matrix back to C++ vector
   for (int i = 0; i < Pmatrix.size(); i++){
    	for (int j = 0; j < Pmatrix.size(); j++){

    		Pmatrix[i][j] = B[count];
    		count++;
    	}
   }

	return Pmatrix;
}
/*
 * Pre-Order traversal works
 *Calculates the JC Matrix
 *
 */
void EvoSim(vector< vector<double> >& rmatrix, Tree * tree, string Ancestor, vector<double>& basefreq){

	double brlength = 0.0;
    vector< vector<double> > QMatrix(4, vector<double>(4, 0.0));
	vector< vector<double> > PMatrix(4, vector<double>(4, 0.0));
	//Pre-Order Traverse the tree
	for (int k = (tree->getNodeCount() - 2); k >= 0; k--){


	    brlength = tree->getNode(k)->getBL();
	    QMatrix = CalQ(rmatrix, basefreq);
	    PMatrix = PCalq(QMatrix, brlength);
	    SeqSim(Ancestor, PMatrix);
	    //If its a tip print the name and the sequence
	    if (tree->getNode(k)->getName() != ""){
	    	string tname = tree->getNode(k)->getName();
	    	//cout << ">" << tname << "\n" << Ancestor << endl;
	    }

	}

}

void randDNA(int len, string& Ancestor,vector<double>& basefreq) {

	string hold = "";
    for (int i = 0; i < len; i++) {
        int RandNumb = 0;
        random_device rand_dev;
        mt19937 generator(rand_dev());
        std::uniform_int_distribution<int>  distr(0, 100);
    	hold = to_string(distr(generator));
    	RandNumb = stod(hold);
        if (RandNumb < 25){
        	Ancestor += "A";
        }else if(RandNumb > 25 && RandNumb < 50){

        	Ancestor += "T";
        }else if(RandNumb > 50 && RandNumb < 75){

        	Ancestor += "C";
        }else{

        	Ancestor += "G";
        }
    }
    cout << Ancestor << endl;
}

//Take in a tree
void SEQGEN::TakeInTree(vector< vector<double> >& rmatrix, Tree * tree, int length, vector<double>& basefreq){

	string Ancestor = "";
	string branch = "";
    vector< vector<double> > Matrix(4, vector<double>(4, 1.0));
    randDNA(length, Ancestor, basefreq);
    EvoSim(rmatrix, tree, Ancestor, basefreq);
}

SEQGEN::SEQGEN() {
	// TODO Auto-generated constructor stub

}

SEQGEN::~SEQGEN() {
	// TODO Auto-generated destructor stub
}
