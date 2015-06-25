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


void SeqSim(string& Ancestor, vector< vector<double> > Matrix){

	//cout << "In Seq Sim" << endl;
    const int range_from  = 0;
    const int range_to    = 100000;
    //Odds of it being any of these states
    float RandNumb = 0.0;
    float ChanceA = 0.0;
    float ChanceT = 0.0;
    float ChanceC = 0.0;
    float ChanceG = 0.0;
    int row = 0;
    string hold = "";
    random_device rand_dev;
    mt19937 generator(rand_dev());
    std::uniform_int_distribution<int>  distr(range_from, range_to);
    string newstring = "";
    for (int i = 0; i < Ancestor.size(); i++){
        float RandNumb = 0.0;
        float ChanceA = 0.0;
        float ChanceT = 0.0;
        float ChanceC = 0.0;
        float ChanceG = 0.0;
        random_device rand_dev;
        mt19937 generator(rand_dev());
        std::uniform_int_distribution<int>  distr(range_from, range_to);
    	hold = to_string(distr(generator));
    	RandNumb = stod(hold);
    	RandNumb /= 100000;
    	if (Ancestor[i] == 'A'){
    		row = 0;
    		ChanceA = abs(Matrix[row][row]);
    		ChanceT = Matrix[row][1] /abs(Matrix[row][row]);
    		ChanceC = (Matrix[row][2] /abs(Matrix[row][row]) + ChanceT);
    		ChanceG = (Matrix[row][3] /abs(Matrix[row][row]) + ChanceC);
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
    		row = 1;
    		ChanceT = abs(Matrix[row][row]);
    		ChanceA = Matrix[row][0] / ChanceT;
    		ChanceC = (Matrix[row][2] / ChanceT + ChanceA);
    		ChanceG = (Matrix[row][3] / ChanceT + ChanceC);
    		if(RandNumb < ChanceA){
    			newstring += "A";
    		}else if(RandNumb > ChanceA & RandNumb < ChanceC){
    			newstring += "C";
    		}else if(RandNumb > ChanceC && RandNumb < ChanceG){

    			newstring += "G";
    		}else{

    			newstring += "T";
    		}
    	}else if (Ancestor[i] == 'C'){
    		row = 2;
    		ChanceC = abs(Matrix[row][row]);
    		ChanceA = Matrix[row][0] / ChanceC;
    		ChanceT = (Matrix[row][1] /ChanceC + ChanceA);
    		ChanceG = (Matrix[row][3] /ChanceC + ChanceT);
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
    		row = 3;
    		ChanceG = abs(Matrix[row][row]);
    		ChanceA = Matrix[row][0] / ChanceG;
    		ChanceT = (Matrix[row][1] / ChanceG + ChanceA);
    		ChanceC = (Matrix[row][2] / ChanceG + ChanceC);
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
 * Calculate the Q Matrix
 *
 *
 */
vector< vector<double> > CalQ(vector< vector<double> > rmatrix, vector<double> basefreq){

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

			//cout << bigpi[i][j] << "\t";
    	}
    	//cout << endl;
    }
    for (int i = 0; i < rmatrix.size(); i++){
    	for (int j = 0; j < rmatrix.size(); j++){

    		if (i != j){

    			bigpi[i][j] /= tscale;

    		}else{
    			bigpi[i][j] = 0.0;
    		}

			//cout << bigpi[i][j] << "\t";
    	}
    	//cout << endl;
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
    for (int i = 0; i < rmatrix.size(); i++){
    	for (int j = 0; j < rmatrix.size(); j++){

			//cout << bigpi[i][j] << "\t";
    	}
    	//cout << endl;
    }

	return bigpi;
}
/*Calculate the P Matrix
 */
vector< vector<double> > PCalq(vector< vector<double> > QMatrix, float br){

	vector<double> PMatrix;
	vector< vector<double> > NewP(4, vector<double>(4, 0.0));
	mat A = randn<mat>(4,4);
	mat B;
   for (int i = 0; i < QMatrix.size(); i++){
    	for (int j = 0; j < QMatrix.size(); j++){

    		if (i == j){

    			//cout << QMatrix[i][j] << " ";

    		}else{
    			//cout << QMatrix[i][j] << " ";
    		}
    		PMatrix.push_back(QMatrix[i][j]);
    	}
    	//cout << endl;
    }
   //Q * t
   for (int i = 0; i < PMatrix.size(); i++){

	   A[i] = (PMatrix[i] * br);
   }
   //e^
   B = expmat(A);
   //cout << B << endl;
   int count = 0;
   for (int i = 0; i < NewP.size(); i++){
    	for (int j = 0; j < NewP.size(); j++){

    		NewP[i][j] = B[count];
    		count++;
    	}
   }

	return NewP;
}
/*
 * Pre-Order traversal works
 *Calculates the JC Matrix
 *
 */
void EvoSim(vector< vector<double> >& rmatrix, Tree * tree, string Ancestor){

	Node nodes;
	double brlength = 0.0;
	string branch = "";
    vector<double> basefreq(4, 0.0);
    vector< vector<double> > QMatrix;
	vector< vector<double> > PMatrix(4, vector<double>(4, 0.0));
    basefreq[0] = 0.25;
    basefreq[1] = 0.25;
    basefreq[2] = 0.25;
    basefreq[3] = 0.25;
	//Pre-Order Traverse the tree
	for (int k = (tree->getNodeCount() - 2); k >= 0; k--){


		//cout << ">node" << to_string(k) << "\n" << Ancestor << endl;
	    brlength = tree->getNode(k)->getBL();
	    QMatrix = CalQ(rmatrix, basefreq);
	    PMatrix = PCalq(QMatrix, brlength);
	    SeqSim(Ancestor, PMatrix);
	    //If their is a name print the name and the sequence
	    if (tree->getNode(k)->getName() != ""){
	    	string tname = tree->getNode(k)->getName();
	    	branch = to_string(brlength);
	    	cout << ">" << tname << "\n" << Ancestor << endl;
	    }

	}

}

void randDNA(int len, string& Ancestor) {
	char DNA[] = "ATCG";

    for (int i = 0; i < len; i++) {
        Ancestor += DNA[rand() % (5 - 1)];
    }
}

//Take in a tree
void SEQGEN::TakeInTree(Tree * tree, int length){

	string Ancestor = "";
	string branch = "";
    vector< vector<double> > Matrix(4, vector<double>(4, 1.0));
    vector< vector<double> > QMatrix;
	vector< vector<double> > PMatrix(4, vector<double>(4, 0.0));


    /*Test Stuff
     * rmatrix is the rate matrix
     * Gets filled with .333 and -1.0
     *To be given in the future
     */
    vector< vector<double> > rmatrix(4, vector<double>(4, 1.0));
    for (int i = 0; i < rmatrix.size(); i++){
    	for (int j = 0; j < rmatrix.size(); j++){

    		if (i == j){//Fill Diagnol

    			rmatrix[i][j] *= -1.0;

    		}else{
    			rmatrix[i][j] *= 0.33333;
    		}
    	}
    }
    randDNA(length, Ancestor);
    EvoSim(rmatrix, tree, Ancestor);
}

SEQGEN::SEQGEN() {
	// TODO Auto-generated constructor stub

}

SEQGEN::~SEQGEN() {
	// TODO Auto-generated destructor stub
}
/*
for(int i=0;i<tree->getExternalNodeCount();i++){
    string tname = tree->getExternalNode(i)->getName();
    brlength = tree->getExternalNode(i)->getBL();
    //JC = ((1.0/4.0)-((1.0/4.0)*exp(-4.0*((brlength))/3.0)));
    //JC = 1/4.0 * (1 - (exp(-4.0*(brlength))));
    branch = to_string(brlength);
    	cout << ">" << tname << "_" << branch << "\n" << Ancestor << endl;
}*/
/*
for(int i=0;i<tree->getNodeCount();i++){
    brlength = tree->getNode(i)->getBL();
    string tname = tree->getNode(i)->getName();
    branch = to_string(brlength);
    cout << branch << "\t" << tname << endl;
}
*/
/*
 * Calculates the rate matrix under JC

void JC(vector< vector<double> >& Matrix, double& brlength){

    vector< vector<double> > TempMatrix(4, vector<double>(4, 1.0));
	//cout << "BR Length: " << brlength << endl;
    for (int i = 0; i < Matrix.size(); i++){
    	for (int j = 0; j < Matrix.size(); j++){

    		if (i == j){//Fill Diagnol

    			TempMatrix[i][j] *= 0.25 + (0.75 * exp((-4.0*brlength) / 3.0));


    		}else{
    			TempMatrix[i][j] *= 0.25 - (0.25 * exp((-4.0*brlength) / 3.0));
    		}
			//cout << TempMatrix[i][j] << "\t";
    	}
    	//cout << endl;
    }
    Matrix = TempMatrix;

}*/
