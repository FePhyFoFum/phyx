#include <string>
#include <vector>
//#include <math.h>
#include <cmath>
#include <armadillo>

using namespace arma;

#include "tree_utils.h"
#include "tree.h"
#include "cont_models.h"
#include "constants.h" // for PI and E


void calc_vcv(Tree * tree, mat & vcv) {
    int numlvs = tree->getExternalNodeCount();
    vcv = mat(numlvs,numlvs);
    int count = 0;
    for (int i=0; i < numlvs; i++) {
        int count2 = 0;
        for (int j=0; j < numlvs; j++) {
            if (i != j) {
                Node * a = getMRCA_forVCV(tree->getExternalNode(i),tree->getExternalNode(j));
                double length = get_length_to_root(a);
                vcv(count,count2) = length;
            } else {
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
Node * getMRCA_forVCV(Node * curn1,Node * curn2) {
    Node * mrca = NULL;
    //get path to root for first node
    vector<Node *> path1;
    Node * parent = curn1;
    path1.push_back(parent);
    while (parent != NULL) {
        path1.push_back(parent);
        if (parent->getParent() != NULL) {
            parent = parent->getParent();
        } else {
            break;
        }
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

Node * getMRCAFromPath_forVCV(vector<Node *> * path1,Node * curn2) {
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
double norm_pdf_multivariate(rowvec & x, rowvec & mu, mat & sigma) {
    unsigned int size = x.n_cols;
    if (size == mu.n_cols && sigma.n_rows == size && sigma.n_cols == size) {
        double DET = det(sigma);
        if (DET == 0) {
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
    } else {
        cerr << "The dimensions of the input don't match" << endl;
        exit(0);
    }
}

double norm_log_pdf_multivariate(rowvec & x, rowvec & mu, mat & sigma) {
    unsigned int size = x.n_cols;
    if (size == mu.n_cols && sigma.n_rows == size && sigma.n_cols == size) {
        double DET;
        double sign;
        log_det(DET,sign,sigma);
        if (DET == 0) {
            cerr << "The covariance matrix can't be singular" << endl;
            exit(0);
        }
        mat U; vec s; mat V;
        svd(U, s, V, sigma,"dc");
        mat diagD (s.size(),s.size());
        diagD.zeros();
        diagD.diag() = 1./s;
        mat invC2 = V*diagD*trans(U);
        rowvec ancA = x - mu;
        double final = -.5 * (dot(invC2*trans(ancA),ancA))-0.5 * DET -0.5 * (size * log(2*PI));
        return final;
    } else {
        cerr << "The dimensions of the input don't match" << endl;
        exit(0);
    }
}

/**
 * assumes that the characters are in get_cont_char and that the 
 * results will be in assocDoubleVector as val and valse
 */
void calc_square_change_anc_states(Tree * tree, int index) {
    int df = 0;
    int count = 0;
    map<Node *,int> nodenum;
    for (int i=0; i < tree->getInternalNodeCount(); i++) {
        nodenum[tree->getInternalNode(i)] = count;
        count += 1;
        df += 1;
        (*tree->getInternalNode(i)->getDoubleVector("val"))[index] = 0.0;
    }
    df -= 1;
    mat fullMcp(df+1,df+1);
    vec fullVcp(df+1);
    fullMcp.fill(0.0);
    fullVcp.fill(0.0);
    calc_postorder_square_change(tree->getRoot(),nodenum,&fullMcp,&fullVcp,index);
    mat b = chol(fullMcp);
    vec mle;
    mat x = solve(trimatl(b.t())*b,fullVcp);
    count = 0;
    for (int i=0; i < tree->getInternalNodeCount(); i++) {
        (*tree->getInternalNode(i)->getDoubleVector("val"))[index] = x(nodenum[tree->getInternalNode(i)],0);
        count += 1;
    }
}

void calc_postorder_square_change(Node * node,map<Node *,int> & nodenum,
    mat * fullMcp, mat * fullVcp, int index) {
    for (int i=0; i < node->getChildCount(); i++) {
        calc_postorder_square_change(node->getChild(i),nodenum,fullMcp,fullVcp,index);    
    }
    if (node->getChildCount() > 0) {
        int nni = nodenum[node];
        for (int j=0; j < node->getChildCount(); j++) {
            double tbl = 2./node->getChild(j)->getBL();
            (*fullMcp)(nni,nni) += tbl;
            if (node->getChild(j)->getChildCount() == 0) {
                (*fullVcp)[nni] += (*node->getChild(j)->getDoubleVector("val"))[index] * tbl;
            } else {
                int nnj = nodenum[node->getChild(j)];
                (*fullMcp)(nni,nnj) -= tbl;
                (*fullMcp)(nnj,nni) -= tbl;
                (*fullMcp)(nnj,nnj) += tbl;
            }
        }
    }
}

double calc_bm_node_postorder(Node * node, int nch, double sigma){
    double node_like = 0.;
    for (int i=0;i<node->getChildCount();i++){
        if(node->getChild(i)->isInternal()){
           node_like += calc_bm_node_postorder(node->getChild(i),nch,sigma);
        }
    }
    if (node->isInternal()){
        double ch1 = (*node->getChild(0)->getDoubleVector("val"))[nch];
        double ch2 = (*node->getChild(1)->getDoubleVector("val"))[nch];
        double ch = ch1 - ch2;
        double bl1 = node->getChild(0)->getBL(); 
        double bl2 = node->getChild(1)->getBL();
        double bl = bl1 + bl2;
        double cur_like = ((-0.5)* ((log(2*M_PI*sigma))+(log(bl))+(pow(ch,2)/(sigma*bl))));
        node_like += cur_like;
        if (node->isRoot() == false){
            (*node->getDoubleVector("val"))[nch] = ((bl2*ch1)+(bl1*ch2))/(bl);
            node->setBL(node->getBL()+((bl1*bl2)/(bl1+bl2)));
        }
    }
    return node_like;
}

double calc_bm_prune(Tree * tr, double sigma){
    int nchar = (*tr->getRoot()->getDoubleVector("val")).size();
    double tlike = 0;
    map<Node *, double> oldlen;
    for (int i=0;i<tr->getNodeCount();i++){oldlen[tr->getNode(i)] = tr->getNode(i)->getBL();}
    for (int i=0;i<nchar;i++){
        for (int j=0;j<tr->getNodeCount();j++){tr->getNode(j)->setBL(oldlen[tr->getNode(j)]);}
        tlike += calc_bm_node_postorder(tr->getRoot(),i,sigma);
    }
    return tlike;
}


