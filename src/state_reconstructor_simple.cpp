
#include <string>
#include <vector>
#include <map>

using namespace std;

#include "state_reconstructor_simple.h"
#include "tree.h"
#include "rate_model.h"
#include "node.h"
#include "utils.h"
#include "vector_node_object.h"
#include "sequence.h"

#define verbose false

#define MINBL 0.000000001

/**
 * this can take the multiple states
 */

StateReconstructorSimple::StateReconstructorSimple(RateModel & _rm, int num_sites):tree(NULL),
    nstates(_rm.nstates),nsites(num_sites),rm(_rm),dc("dist_conditionals"),v_storage(_rm.nstates,0),
    v1(_rm.nstates,0),v2(_rm.nstates,0) {}
    /*
     * initialize each node with segments
     */
void StateReconstructorSimple::set_tree(Tree * tr) {
    tree = tr;
    if (verbose) {
        cout << "initializing nodes..." << endl;
    }
    for (int i=0; i < tree->getNodeCount(); i++) {
        if (tree->getNode(i)->getBL()<MINBL) {
            tree->getNode(i)->setBL(MINBL * 100);
    }
    //VectorNodeObject<double> * dcs = new VectorNodeObject<double>(nstates);
    vector<vector<double> > tv;
    for (int j=0; j < nsites; j++) {
        vector<double> tn(nstates);
        tv.push_back(tn);
    }
    conditionals_map[tree->getNode(i)] = tv;
    //tree->getNode(i)->assocObject(dc,*dcs);
    //delete dcs;
    }
}

bool StateReconstructorSimple::set_tip_conditionals(vector<Sequence> & distrib_data, int site) {
    bool allsame = true;
    string testsame = distrib_data[0].get_sequence();
    for (unsigned int i=0; i < distrib_data.size(); i++) {
        Sequence seq = distrib_data[i];
        Node * nd = tree->getExternalNode(seq.get_id());
        if (verbose) {
            cout << nd->getName() << " ";
        }
        for (int j=0; j < nstates; j++) {
            if (seq.get_sequence().at(j) == '1') {
                //(((VectorNodeObject<double>*) nd->getObject(dc)))->at(j) = 1.0;
                conditionals_map[nd][site][j] = 1.0;
            } else {
                //(((VectorNodeObject<double>*) nd->getObject(dc)))->at(j) = 0.0;
                conditionals_map[nd][site][j] = 0.0;
            }
            if (verbose) {
                cout << seq.get_sequence().at(j);
            }
        }
        if (verbose) {
            cout << endl;
        }
        if (testsame != seq.get_sequence()) {
            allsame = false;
        }
    }
    if (allsame == true && verbose == true) {
        cerr << "all the tips have the same characters" << endl;
    }
    return allsame;
}

//VectorNodeObject<double> 
void StateReconstructorSimple::conditionals(vector<double> * v, Node & node,int site) {
//    VectorNodeObject<double> distconds = *((VectorNodeObject<double>*) node.getObject(dc));
    vector<double> * tdistconds = &conditionals_map[&node][site];

    //VectorNodeObject<double> * v = new VectorNodeObject<double> (nstates, 0);
    if (map_ps.count(node.getBL()) == 0) {
    mat tp(nstates,nstates);
    rm.setup_fortran_P(tp,node.getBL(),false);
    map_ps[node.getBL()] = tp;
    }
    mat * p = &map_ps[node.getBL()];
    for (int j=0; j < nstates; j++) {
        v_storage[j] = 0;
        for ( int k=0; k < nstates; k++) {
            //v_storage[j] += (distconds.at(k)*(*p)(j,k));
            v_storage[j] += tdistconds->at(k)*(*p)(j, k);
        }
    }
    for (unsigned int j=0; j < tdistconds->size(); j++) {
        //distconds[j] = v_storage[j];
        //tdistconds->at(j) = v_storage[j];
        (*v)[j] = v_storage[j];
    }
//    return distconds;
}

//model 2a
void StateReconstructorSimple::conditionals2(vector<double> * v, Node & node,int site) {
    vector<double> * tdistconds = &conditionals_map[&node][site];
    if (map_ps0.count(node.getBL()) == 0) {
        mat tp0(nstates,nstates);
        mat tp1(nstates,nstates);
        mat tp2(nstates,nstates);
        rm.setup_fortran_P_whichQ(0,tp0,node.getBL());
        rm.setup_fortran_P_whichQ(1,tp1,node.getBL());
        rm.setup_fortran_P_whichQ(2,tp2,node.getBL());

        map_ps0[node.getBL()] = tp0;
        map_ps1[node.getBL()] = tp1;
        map_ps2[node.getBL()] = tp2;
    }
    mat * p0 = &map_ps0[node.getBL()];
    mat * p1 = &map_ps1[node.getBL()];
    mat * p2 = &map_ps2[node.getBL()];

    for (int j=0; j < nstates; j++) {
        v_storage[j] = 0;
        for (int k=0; k < nstates; k++) {
            //v_storage[j] += (distconds.at(k)*(*p)(j, k));
            v_storage[j] += tdistconds->at(k)*(*p0)(j, k) * pp0;
            v_storage[j] += tdistconds->at(k)*(*p1)(j, k) * pp1;
            v_storage[j] += tdistconds->at(k)*(*p2)(j, k) * pp2;
        }
    }
    for (unsigned int j=0; j < tdistconds->size(); j++) {
        (*v)[j] = v_storage[j];
    }
}

void StateReconstructorSimple::ancdist_conditional_lh(Node & node, int site) {
//    VectorNodeObject<double> distconds(nstates, 0);
    vector<double> * distconds = &conditionals_map[&node][site];
    if (node.isExternal() == false) {//is not a tip
        Node * c1 = node.getChild(0);
        Node * c2 = node.getChild(1);
        ancdist_conditional_lh(*c1, site);
        ancdist_conditional_lh(*c2, site);
        //VectorNodeObject<double> v1;
        //VectorNodeObject<double> v2;
        //v1 =conditionals(*c1,site);
        //v2 =conditionals(*c2,site);
        if (rm.selection_model == 0) {
            conditionals(&v1,*c1, site);
            conditionals(&v2,*c2, site);
        } else if (rm.selection_model == 2) {
            conditionals2(&v1,*c1, site);
            conditionals2(&v2,*c2, site);
        }
        //vector<Superdouble> * v1 = &conditionals_map[c1][site];
        //vector<Superdouble> * v2 = &conditionals_map[c2][site];
        for (int i=0; i < nstates; i++) {
    //        distconds.at(i)= v1[i] * v2[i];
            distconds->at(i) = v1[i] * v2[i];
        }
        if (node.isRoot()) {
            //with equal freq!
            for (int i=0; i < nstates; i++) {
                distconds->at(i) = distconds->at(i) * (1./nstates);
            }
        }
    }
//    } else {
//    distconds = *((VectorNodeObject<double>*)node.getObject(dc));
//    }
//    for (unsigned int i=0; i < distconds.size(); i++) {
//    ((VectorNodeObject<double>*)node.getObject(dc))->at(i) = distconds.at(i);
//    }
}

double StateReconstructorSimple::eval_likelihood(int site) {
   ancdist_conditional_lh(*tree->getRoot(), site);
//    return (-log(sum(*
//          (VectorNodeObject<double>*) tree->getRoot()->getObject(dc))));
    return (-log(sum(conditionals_map[tree->getRoot()][site])));
    //return double(-(calculate_vector_Superdouble_sum(*(VectorNodeObject<double>*) tree->getRoot()->getObject(dc))).getLn());
}

void StateReconstructorSimple::clear_map_ps() {
    if (rm.selection_model == 0) {
        map_ps.clear();
    } else if (rm.selection_model == 1) {
        
        // something supposed to go here?
        
    } else if (rm.selection_model == 2) {
        map_ps0.clear();
        map_ps1.clear();
        map_ps2.clear();
    }
}


StateReconstructorSimple::~StateReconstructorSimple() {

}

