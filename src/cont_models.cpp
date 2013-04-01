#include <string>
#include <vector>
#include <math.h>
#include "tree_utils.h"
#include "tree.h"
#include "cont_models.h"

#include <armadillo>
using namespace arma;

#define PI 3.1415926535897932384626433832795
#define E 2.718281828459045


void calc_vcv(Tree * tree, mat & vcv){
    int numlvs = tree->getExternalNodeCount();
    vcv = mat(numlvs,numlvs);
    int count = 0;
    for (int i=0;i<numlvs;i++){
        int count2 = 0;
        for (int j=0;j<numlvs;j++){
            if (i != j){
                Node * a = getMRCA_forVCV(tree->getExternalNode(i),tree->getExternalNode(j));
                double length = get_length_to_root(a);
                vcv(count,count2) = length;
            }else{
                double length = get_length_to_root(tree->getExternalNode(i));
                vcv(count,count2) = length;
            }
            count2 += 1;
        }
        count += 1;
    }
}

/*
 * get the MRCA
 * this calculates the typical algorithm for MRCA
 * can be a little slow, so probably best to use
 * getMRCAFromPath_forVCV
 */
Node * getMRCA_forVCV(Node * curn1,Node * curn2){
    Node * mrca = NULL;
    //get path to root for first node
    vector<Node *> path1;
    Node * parent = curn1;
    path1.push_back(parent);
    while (parent != NULL) {
        path1.push_back(parent);
        if (parent->getParent() != NULL)
            parent = parent->getParent();
        else
            break;
    }
    //find first match between this node and the first one
    parent = curn2;
    bool x = true;
    while (x == true) {
        int psize = path1.size();
        for (int i = 0; i < psize; i++) {
            if (parent == path1.at(i)) {
                mrca = parent;
                x = false;
                break;
            }
        }
        parent = parent->getParent();
    }
    return mrca;
}

/*
 * this saves a great deal of time as it takes the already
 * obtained path and finds the match with the second node
 */

Node * getMRCAFromPath_forVCV(vector<Node *> * path1,Node * curn2){
    Node * mrca = NULL;
    Node * parent = curn2;
    bool x = true;
    while (x == true) {
        int psize = path1->size();
        for (int i = 0; i < psize; i++) {
            if (parent == path1->at(i)) {
                mrca = parent;
                x = false;
                break;
            }
        }
        parent = parent->getParent();
    }
    return mrca;
}


/*
   calculate the maximum likelihood for the ancestral state and the rate
   calc_bm_like only calculates the rate and solves for the anc state
#http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Non-degenerate_case
x should be a vector
mu should be a vector
sigma should be vcv with sigma already applied
 */
double norm_pdf_multivariate(rowvec & x, rowvec & mu, mat & sigma){
    int size = x.n_cols;
    if (size == mu.n_cols && sigma.n_rows == size && sigma.n_cols == size){
        double DET = det(sigma);
        if (DET == 0){
            cerr << "The covariance matrix can't be singular" << endl;
            exit(0);
        }
        double norm_const = 1.0/ ( pow((2*PI),double(size)/2) * pow(DET,1.0/2) );
        mat x_mu = (mat) x - mu;
        mat tm = (x_mu * inv(sigma) * trans(x_mu));
        double big1 = tm(0,0);
        double result = pow(E, -0.5 * big1);
        double final = norm_const * result;
        return final;
    }else{
        cerr << "The dimensions of the input don't match" << endl;
        exit(0);
    }
}

double norm_log_pdf_multivariate(rowvec & x, rowvec & mu, mat & sigma){
    int size = x.n_cols;
    if (size == mu.n_cols && sigma.n_rows == size && sigma.n_cols == size){
        double DET;
        double sign;
        log_det(DET,sign,sigma);
        if (DET == 0){
            cerr << "The covariance matrix can't be singular" << endl;
            exit(0);
        }
        mat U; vec s; mat V;
        svd(U, s, V, sigma);//, method = "dc") 
        mat diagD (s.size(),s.size());
        diagD.zeros();
        diagD.diag() = 1./s;
        mat invC2 = V*diagD*trans(U);
        rowvec ancA = x - mu;
        double final = -.5 * (dot(invC2*trans(ancA),ancA))-0.5 * DET -0.5 * (size * log(2*PI));
        return final;
    }else{
        cerr << "The dimensions of the input don't match" << endl;
        exit(0);
    }
}

