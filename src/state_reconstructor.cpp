
#include <string>
#include <vector>
#include <map>

using namespace std;

#include "state_reconstructor.h"
#include "tree.h"
#include "rate_model.h"
#include "node.h"
#include "utils.h"
#include "vector_node_object.h"
#include "sequence.h"
#include "superdouble.h"

#define verbose true

#define MINBL 0.000000001

StateReconstructor::StateReconstructor(RateModel & _rm,vector<RateModel> &_vrm):tree(NULL),
    use_periods(false),nstates(_rm.nstates),rm(_rm),rm_periods(_vrm),dc("dist_conditionals"),
    andc("anc_dist_conditionals"),store_p_matrices(false),use_stored_matrices(false),revB("revB"),
    rev(false),rev_exp_number("rev_exp_number"),rev_exp_time("rev_exp_time"),
    stochastic(false),stored_EN_matrices(map<Superdouble, mat >()),
    stored_ER_matrices(map<Superdouble, mat >()),sp_alphas("sp_alphas"),alphas("alphas") {}

/**
 * need to do this before you do the set tree
 */
void StateReconstructor::set_periods(vector<double> & ps,vector<RateModel> & rms) {
    use_periods = true;
    periods = ps;
    rm_periods = rms;
}

/*
 * initialize each node with segments
 */
void StateReconstructor::set_tree(Tree * tr) {
    tree = tr;
    if (verbose) {
        cout << "initializing nodes..." << endl;
    }
    for (int i=0; i < tree->getNodeCount(); i++) {
        if (tree->getNode(i)->getBL()<MINBL) {
            tree->getNode(i)->setBL(MINBL * 100);
        }
        if (use_periods == true) {
            tree->getNode(i)->initSegVector();
        }
        VectorNodeObject<Superdouble> * dcs = new VectorNodeObject<Superdouble>(nstates);
        tree->getNode(i)->assocObject(dc,*dcs);
        delete dcs;
    }
    /*
     * initialize the actual branch segments for each node
     */
    if (use_periods == true) {
        cout << "initializing branch segments..." << endl;
        tree->setHeightFromTipToNodes();
        for (int i=0; i < tree->getNodeCount(); i++) {
            if (tree->getNode(i)->hasParent()) {
                vector<double> pers(periods);
                double anc = tree->getNode(i)->getParent()->getHeight();
                double des = tree->getNode(i)->getHeight();
                double t = des;
                if (pers.size() > 0) {
                    for (unsigned int j=0; j < pers.size(); j++) {
                        double s = 0;
                        if (pers.size() == 1)
                            s = pers[0];
                        for (unsigned int k=0; k < j+1; k++) {
                            s += pers[k];
                        }
                        if (t < s) {
                            double duration = min(s-t,anc-t);
                            if (duration > 0) {
                                BranchSegment tseg = BranchSegment(duration,j);
                                tree->getNode(i)->getSegVector()->push_back(tseg);
                            }
                            t += duration; // TODO: make sure that this is all working
                        }
                        if (t > anc || pers[j] > t) {
                            break;
                        }
                    }
                } else {
                    BranchSegment tseg = BranchSegment(tree->getNode(i)->getBL(),0);
                    tree->getNode(i)->getSegVector()->push_back(tseg);
                }
            }
        }
    }
}

/**
 * this will setup the distconds and ancdistconds for each segment
 */
void StateReconstructor::set_periods_model() {
    for (int i=0; i < tree->getNodeCount(); i++) {
        vector<BranchSegment> * tsegs = tree->getNode(i)->getSegVector();
        for (unsigned int j=0; j < tsegs->size(); j++) {
            tsegs->at(j).setModel(&rm_periods[tsegs->at(j).getPeriod()]);
            vector<Superdouble> * distconds = new vector<Superdouble> (nstates,0);
            tsegs->at(j).distconds = distconds;
            vector<Superdouble> * ancdistconds = new vector<Superdouble> (nstates,0);
            tsegs->at(j).ancdistconds = ancdistconds;
        }
    }
    vector<Superdouble> * distconds = new vector<Superdouble> (nstates,0);
    tree->getRoot()->assocDoubleVector(dc,*distconds);
    delete distconds;
    vector<Superdouble> * ancdistconds = new vector<Superdouble> (nstates,0);
    tree->getRoot()->assocDoubleVector(andc, *ancdistconds);
    delete ancdistconds;
}

bool StateReconstructor::set_tip_conditionals(vector<Sequence> & data) {
    bool allsame = true;
    string testsame = data[0].get_sequence();
    for (unsigned int i=0; i < data.size(); i++) {
        Sequence seq = data[i];
        Node * nd = tree->getExternalNode(seq.get_id());
        if (verbose) {
            cout << nd->getName() << " ";
            }
        if (use_periods == false) {
            for (int j=0; j < nstates; j++) {
            if (seq.get_sequence().at(j) == '1') {
                (((VectorNodeObject<Superdouble>*) nd->getObject(dc)))->at(j) = 1.0;
            } else {
                (((VectorNodeObject<Superdouble>*) nd->getObject(dc)))->at(j) = 0.0;
                    }
            if (verbose) {
                cout << seq.get_sequence().at(j);
                    }
            }
        } else {
            vector<BranchSegment> * tsegs = nd->getSegVector();
            for (int j=0; j < nstates; j++) {
                if (seq.get_sequence().at(j) == '1') {
                    tsegs->at(0).distconds->at(j) = 1.0;
                } else {
                    tsegs->at(0).distconds->at(j) = 0.0;
                        }
                if (verbose) {
                    cout << seq.get_sequence().at(j);
                }
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

bool StateReconstructor::set_tip_conditionals_already_given(vector<Sequence> & data) {
    bool allsame = false;
    for (unsigned int i=0; i < data.size(); i++) {
    Sequence seq = data[i];
    Node * nd = tree->getExternalNode(seq.get_id());
    vector<string> searchtokens;
    tokenize(seq.get_sequence(), searchtokens, ",");
    for (unsigned int j=0; j < searchtokens.size(); j++) {
        trim_spaces(searchtokens[j]);
    }
    if (verbose) {
        cout << nd->getName() << " ";
    }
    if (use_periods == false) {
        for (int j=0; j < nstates; j++) {
            (((VectorNodeObject<Superdouble>*) nd->getObject(dc)))->at(j) = atof(searchtokens[j].c_str());
            if (verbose) {
                cout << searchtokens[j];
            }
        }
    } else {
        vector<BranchSegment> * tsegs = nd->getSegVector();
        for (int j=0; j < nstates; j++) {
        tsegs->at(0).distconds->at(j) = atof(searchtokens[j].c_str());
        if (verbose) {
            cout << searchtokens[j];
                }
        }
    }    
    if (verbose) {
        cout << endl;
        }
    }
    return allsame;
}


VectorNodeObject<Superdouble> StateReconstructor::conditionals(Node & node) {
    VectorNodeObject<Superdouble> distconds = *((VectorNodeObject<Superdouble>*) node.getObject(dc));
    VectorNodeObject<Superdouble> * v = new VectorNodeObject<Superdouble> (nstates, 0);
    cx_mat p;
    if (use_stored_matrices == false) {
        p= rm.setup_P(node.getBL(), store_p_matrices);
    } else {
        p = rm.stored_p_matrices[node.getBL()];
    }
    for (int j=0; j < nstates; j++) {
        for (int k=0; k < nstates; k++) {
            v->at(j) += (distconds.at(k)*real(p(j, k)));
        }
    }
    for (unsigned int j=0; j < distconds.size(); j++) {
        distconds[j] = v->at(j);
    }
    if (store_p_matrices == true) {
        node.assocObject(sp_alphas,distconds);
        node.assocObject(alphas,distconds);
    }
    delete v;
    return distconds;
}

VectorNodeObject<Superdouble> StateReconstructor::conditionals_periods(Node & node) {
    vector<Superdouble> distconds;
    vector<BranchSegment> * tsegs = node.getSegVector();
    distconds = *tsegs->at(0).distconds;
    for (unsigned int i=0; i < tsegs->size(); i++) {
        for (unsigned int j=0; j < distconds.size(); j++) {
            tsegs->at(i).distconds->at(j) = distconds.at(j);
        }
        RateModel * trm = tsegs->at(i).getModel();
        vector<Superdouble> * v = new vector<Superdouble> (nstates, 0);
        //vector<vector<double > > p;
        cx_mat p;
        if (use_stored_matrices == false) {
            //p= trm->setup_fortran_P(tsegs->at(i).getPeriod(),tsegs->at(i).getDuration(),store_p_matrices);
            p = trm->setup_P(tsegs->at(i).getDuration(), store_p_matrices);
        } else {
            //p = trm->stored_p_matrices[tsegs->at(i).getPeriod()][tsegs->at(i).getDuration()];
            p = trm->stored_p_matrices[tsegs->at(i).getDuration()];
        }
        for (int j=0; j < nstates; j++) {
            for (int k=0; k < nstates; k++) {
            v->at(j) += (distconds.at(k)*real(p(j, k)));
            }
        }

        for (int j=0; j < nstates; j++) {
            distconds[j] = v->at(j);
        }
        if (store_p_matrices == true) {
            tsegs->at(i).seg_sp_alphas = distconds;
        }
        delete v;
    }
    /*
     * if store is true we want to store the conditionals for each node
     * for possible use in ancestral state reconstruction
     */
    if (store_p_matrices == true) {
        tsegs->at(0).alphas = distconds;
    }
    VectorNodeObject<Superdouble> rdistconds(distconds.size());
    for (unsigned int i=0; i < distconds.size(); i++) {
        rdistconds[i] = distconds[i];
    }
    return rdistconds;
}

void StateReconstructor::ancdist_conditional_lh(Node & node) {
    VectorNodeObject<Superdouble> distconds(nstates, 0);
    if (node.isExternal() == false) { // is not a tip
        Node * c1 = node.getChild(0);
        Node * c2 = node.getChild(1);
        ancdist_conditional_lh(*c1);
        ancdist_conditional_lh(*c2);
        VectorNodeObject<Superdouble> v1;
        VectorNodeObject<Superdouble> v2;
        if (use_periods == false) {
            v1 =conditionals(*c1);
            v2 =conditionals(*c2);
        } else {
            v1 = conditionals_periods(*c1);
            v2 = conditionals_periods(*c2);
        }
        for (int i=0; i < nstates; i++) {
            distconds.at(i)= v1[i] * v2[i];
        }
        //if (node.isRoot()) {
            //with equal freq!
        //    for (int i=0; i < nstates; i++) {
    //        distconds.at(i) = distconds.at(i) * (1./nstates);
    //        }
        //}
    } else {
        if (use_periods == false) {
            distconds = *((VectorNodeObject<Superdouble>*)node.getObject(dc));
        } else {
            vector<BranchSegment> * tsegs = node.getSegVector();
            //distconds = *tsegs->at(0).distconds;
            for (unsigned int i=0; i < distconds.size(); i++) {
                distconds[i] = tsegs->at(0).distconds->at(i);
            }
        }
    }
    if (use_periods == false) {
        for (unsigned int i=0; i < distconds.size(); i++) {
            ((VectorNodeObject<Superdouble>*)node.getObject(dc))->at(i) = distconds.at(i);
        }
    } else {
        if (node.hasParent() == true) {
            vector<BranchSegment> * tsegs = node.getSegVector();
            for (unsigned int i=0; i < distconds.size(); i++) {
                tsegs->at(0).distconds->at(i) = distconds.at(i);
            }
        } else {
            for (unsigned int i=0; i < distconds.size(); i++) {
                ((VectorNodeObject<Superdouble>*)node.getObject(dc))->at(i) = distconds.at(i);
            }
        }
    }
}

double StateReconstructor::eval_likelihood() {
    ancdist_conditional_lh(*tree->getRoot());
    //return (-log(calculate_vector_double_sum(*
    //      (VectorNodeObject<Superdouble>*) tree->getRoot()->getObject(dc))));
    return double(-(calculate_vector_Superdouble_sum(*(VectorNodeObject<Superdouble>*) tree->getRoot()->getObject(dc))).getLn());
}

void StateReconstructor::prepare_ancstate_reverse() {
    reverse(tree->getRoot());
}

void StateReconstructor::reverse(Node * node) {
    rev = true;
    VectorNodeObject<Superdouble> * revconds = new VectorNodeObject<Superdouble> (nstates, 0);//need to delete this at some point
    if (node == tree->getRoot()) {
        for (int i=0; i < nstates; i++) {
            revconds->at(i) = 1.0;//prior
        }
        node->assocObject(revB,*revconds);
        delete revconds;
        for (int i = 0; i < node->getChildCount(); i++) {
            reverse(node->getChild(i));
        }
    } else {
        //else if (node.isExternal() == false) {
        //calculate A i
        //sum over all alpha k of sister node of the parent times the priors of the speciations
        //(weights) times B of parent j
        VectorNodeObject<Superdouble> * parrev = ((VectorNodeObject<Superdouble>*)node->getParent()->getObject(revB));
        VectorNodeObject<Superdouble> sisdistconds;
        if (node->getParent()->getChild(0) != node) {
            VectorNodeObject<Superdouble>* talph = ((VectorNodeObject<Superdouble>*) node->getParent()->getChild(0)->getObject(alphas));
            sisdistconds = *talph;
        } else {
            VectorNodeObject<Superdouble>* talph = ((VectorNodeObject<Superdouble>*) node->getParent()->getChild(1)->getObject(alphas));
            sisdistconds = *talph;
        }

        VectorNodeObject<Superdouble> tempA (nstates,0);
        //needs to be the same as ancdist_cond_lh
        for (int i = 0; i < nstates; i++) {
            //root has i, curnode has left, sister of cur has right
            //for (int j = 0; j < nstates; j++) {
            tempA[i] += (sisdistconds.at(i)*parrev->at(i));
            //}
        }
        //now calculate node B
        //VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
        for (int j=0; j < nstates; j++) {
            revconds->at(j) = 0;
        }
        //RateModel * rm = tsegs->at(ts).getModel();
        cx_mat * p = &rm.stored_p_matrices[node->getBL()];
        mat * EN = NULL;
        mat * ER = NULL;
        VectorNodeObject<Superdouble> tempmoveAer(tempA);
        VectorNodeObject<Superdouble> tempmoveAen(tempA);
        if (stochastic == true) {
            //initialize the segment B's
            for (int j=0; j < nstates; j++) {
                tempmoveAer[j] = 0;
                tempmoveAen[j] = 0;
            }
            EN = &stored_EN_matrices[node->getBL()];
            ER = &stored_ER_matrices[node->getBL()];
        }
        for (int j=0; j < nstates; j++) {
            for (int i = 0; i < nstates; i++) {
                revconds->at(j) += tempA[i]*real((*p)(i, j));//tempA needs to change each time
                if (stochastic == true) {
                    tempmoveAer[j] += tempA[i]*(((*ER)(i, j)));
                    tempmoveAen[j] += tempA[i]*(((*EN)(i, j)));
                }
            }
        }
        for (int j=0; j < nstates; j++) {
            tempA[j] = revconds->at(j);
        }
        if (stochastic == true) {
            node->seg_sp_stoch_map_revB_time = tempmoveAer;
            node->seg_sp_stoch_map_revB_number = tempmoveAen;
        }
        node->assocObject(revB,*revconds); // leak
        delete revconds;
        for (int i = 0; i < node->getChildCount(); i++) {
            reverse(node->getChild(i));
        }
    }
}

vector<Superdouble> StateReconstructor::calculate_ancstate_reverse_sd(Node & node) {
    vector<Superdouble> LHOODS (nstates,0);
    if (node.isExternal() == false) {//is not a tip
        VectorNodeObject<Superdouble> * Bs = (VectorNodeObject<Superdouble> *) node.getObject(revB);
        Node * c1 = node.getChild(0);
        Node * c2 = node.getChild(1);
        VectorNodeObject<Superdouble>* v1  = ((VectorNodeObject<Superdouble>*) c1->getObject(alphas));
        VectorNodeObject<Superdouble>* v2 = ((VectorNodeObject<Superdouble>*) c2->getObject(alphas));

        for (int i = 0; i < nstates; i++) {
            //for (int j=0; j < nstates; j++) {
            //  LHOODS[i] += (v1->at(i)*v2->at(j));//*weight);
            //}
            LHOODS[i] = ((v1->at(i)*v2->at(i)) * Bs->at(i));
        }
    }
    return LHOODS;
}

vector<double> StateReconstructor::calculate_ancstate_reverse(Node & node) {
    vector<double> LHOODS (nstates,0);
    if (node.isExternal() == false) {//is not a tip
        VectorNodeObject<Superdouble> * Bs = (VectorNodeObject<Superdouble> *) node.getObject(revB);
        Node * c1 = node.getChild(0);
        Node * c2 = node.getChild(1);
        VectorNodeObject<Superdouble>* v1  = ((VectorNodeObject<Superdouble>*) c1->getObject(alphas));
        VectorNodeObject<Superdouble>* v2 = ((VectorNodeObject<Superdouble>*) c2->getObject(alphas));

        for (int i = 0; i < nstates; i++) {
            //for (int j=0; j < nstates; j++) {
            //  LHOODS[i] += (v1->at(i)*v2->at(j));//*weight);
            //}
            LHOODS[i] = double((v1->at(i)*v2->at(i)) * Bs->at(i));
        }
    }
    return LHOODS;
}

void StateReconstructor::prepare_stochmap_reverse_all_nodes(int from, int to) {
    stochastic = true;
    //calculate and store local expectation matrix for each branch length
    for (int k = 0; k < tree->getNodeCount(); k++) {
        double dur =  tree->getNode(k)->getBL();
        cx_mat eigvec(nstates, nstates); eigvec.fill(0);
        cx_mat eigval(nstates, nstates); eigval.fill(0);
        bool isImag = rm.get_eigenvec_eigenval_from_Q(&eigval, &eigvec);
        mat Ql(nstates, nstates);Ql.fill(0); Ql(from,to) = rm.get_Q()(from,to);
        mat W(nstates, nstates);W.fill(0); W(from,from) = 1;
        cx_mat summed(nstates, nstates); summed.fill(0);
        cx_mat summedR(nstates, nstates); summedR.fill(0);
        for (int i=0; i < nstates; i++) {
            mat Ei(nstates, nstates); Ei.fill(0);Ei(i,i)=1;
            cx_mat Si(nstates, nstates);
            Si = eigvec * Ei * inv(eigvec);
            for (int j=0; j < nstates; j++) {
                cx_double dij = (eigval(i,i)-eigval(j, j)) * dur;
                mat Ej(nstates, nstates); Ej.fill(0); Ej(j, j)=1;
                cx_mat Sj(nstates, nstates);
                Sj = eigvec * Ej * inv(eigvec);
                cx_double Iijt = 0;
                if (abs(dij) > 10) {
                    Iijt = (exp(eigval(i,i)*dur)-exp(eigval(j, j)*dur))/(eigval(i,i)-eigval(j, j));
                } else if (abs(dij) < 10e-20) {
                    Iijt = dur*exp(eigval(j, j)*dur)*(1.+dij/2.+pow(dij,2.)/6.+pow(dij,3.)/24.);
                } else {
                    if (eigval(i,i) == eigval(j, j)) {
                        //WAS Iijt = dur*exp(eigval(j, j)*dur)*expm1(dij)/dij;
                        if (isImag) {
                            Iijt = dur*exp(eigval(j, j)*dur)*(exp(dij)-1.)/dij;
                        } else {
                            Iijt = dur*exp(eigval(j, j)*dur)*(expm1(real(dij)))/dij;
                        }
                    } else {
                        //WAS Iijt = -dur*exp(eigval(i,i)*dur)*expm1(-dij)/dij;
                        if (isImag) {
                            Iijt = -dur*exp(eigval(i,i)*dur)*(exp(-dij)-1.)/dij;
                        } else {
                            Iijt = -dur*exp(eigval(i,i)*dur)*(expm1(real(-dij)))/dij;
                        }
                    }
                }
                summed += (Si  * Ql * Sj * Iijt);
                summedR += (Si * W * Sj * Iijt);
            }
        }
        //cout << isImag << endl;
        //cout << summed << endl;
        //seems like when these are IMAG, there can sometimes be negative with very small values
        stored_EN_matrices[dur] = abs(real(summed)); //(real(summed));
        stored_ER_matrices[dur] = abs(real(summedR)); //(real(summedR));
    }
}

/*
 * only for number of changes
 */
void StateReconstructor::prepare_stochmap_reverse_all_nodes_all_matrices() {
    stochastic = true;
    //calculate and store local expectation matrix for each branch length
    for (int k = 0; k < tree->getNodeCount(); k++) {
        double dur =  tree->getNode(k)->getBL();
        cx_mat eigvec(nstates, nstates);eigvec.fill(0);
        cx_mat eigval(nstates, nstates);eigval.fill(0);
        bool isImag = rm.get_eigenvec_eigenval_from_Q(&eigval, &eigvec);
        mat Ql(nstates, nstates);
        Ql.fill(0);
        for (unsigned int i=0; i < Ql.n_rows; i++) {
            for (unsigned int j=0; j < Ql.n_cols; j++) {
                if (i != j) {
                    Ql(i, j) = rm.get_Q()(i, j);
                }
            }
        }
        mat W(nstates, nstates); W.fill(0); W(1,1) = 1;
        cx_mat summed(nstates, nstates); summed.fill(0);
        cx_mat summedR(nstates, nstates); summedR.fill(0);
        for (int i=0; i < nstates; i++) {
            mat Ei(nstates, nstates);Ei.fill(0);Ei(i,i)=1;
            cx_mat Si(nstates, nstates);
            Si = eigvec * Ei * inv(eigvec);
            for (int j=0; j < nstates; j++) {
                cx_double dij = (eigval(i,i)-eigval(j, j)) * dur;
                mat Ej(nstates, nstates);Ej.fill(0);Ej(j, j)=1;
                cx_mat Sj(nstates, nstates);
                Sj = eigvec * Ej * inv(eigvec);
                cx_double Iijt = 0;
                if (abs(dij) > 10) {
                    Iijt = (exp(eigval(i,i)*dur)-exp(eigval(j, j)*dur))/(eigval(i,i)-eigval(j, j));
                } else if (abs(dij) < 10e-20) {
                    Iijt = dur*exp(eigval(j, j)*dur)*(1.+dij/2.+pow(dij,2.)/6.+pow(dij,3.)/24.);
                } else {
                    if (eigval(i,i) == eigval(j, j)) {
                        //WAS Iijt = dur*exp(eigval(j, j)*dur)*expm1(dij)/dij;
                        if (isImag) {
                            Iijt = dur*exp(eigval(j, j)*dur)*(exp(dij)-1.)/dij;
                        } else {
                            Iijt = dur*exp(eigval(j, j)*dur)*(expm1(real(dij)))/dij;
                        }
                    } else {
                        //WAS Iijt = -dur*exp(eigval(i,i)*dur)*expm1(-dij)/dij;
                        if (isImag) {
                            Iijt = -dur*exp(eigval(i,i)*dur)*(exp(-dij)-1.)/dij;
                        } else {
                            Iijt = -dur*exp(eigval(i,i)*dur)*(expm1(real(-dij)))/dij;
                        }
                    }
                }
                summed += (Si  * Ql * Sj * Iijt);
                summedR += (Si * W * Sj * Iijt);
            }
        }
        //cout << isImag << endl;
        //cout << summed << endl;
        //seems like when these are IMAG, there can sometimes be negative with very small values
        stored_EN_matrices[dur] = abs(real(summed));//(real(summed));
        stored_ER_matrices[dur] = abs(real(summedR));;//(real(summedR));
    }
}


vector<double> StateReconstructor::calculate_reverse_stochmap(Node & node, bool tm) {
    if (node.isExternal()==false) {//is not a tip
    vector<double> totalExp (nstates,0);
    vector<Superdouble> Bs;
    if (tm) {
        Bs = node.seg_sp_stoch_map_revB_time;
    } else {
        Bs =  node.seg_sp_stoch_map_revB_number;
    }
    Node * c1 = node.getChild(0);
    Node * c2 = node.getChild(1);
    VectorNodeObject<Superdouble> * v1 = ((VectorNodeObject<Superdouble>*) c1->getObject(alphas));
    VectorNodeObject<Superdouble> * v2 = ((VectorNodeObject<Superdouble>*) c2->getObject(alphas));
    VectorNodeObject<double> LHOODS (nstates,0);
    for (int i = 0; i < nstates; i++) {
        //for (int j=0; j < nstates; j++) {
        //int ind1 = leftdists[j];
        //int ind2 = rightdists[j];
        //LHOODS[i] += (v1.at(ind1)*v2.at(ind2)*weight);
        //}
        LHOODS[i] = double(v1->at(i) * v2->at(i) * Bs.at(i));
        //cout << v1->at(i) << " " <<  v2->at(i)<< " " << Bs.at(i) << endl;
    }
    for (int i=0; i < nstates; i++) {
        totalExp[i] = LHOODS[i];
    }
    //not sure if this should return a Superdouble or not when doing a bigtree
    return totalExp;
    } else {
    vector<double> totalExp (nstates,0);
    vector<Superdouble> Bs;
    if (tm) {
        Bs = node.seg_sp_stoch_map_revB_time;
    } else {
        Bs =  node.seg_sp_stoch_map_revB_number;
        }
    VectorNodeObject<double> LHOODS (nstates,0);
    VectorNodeObject<Superdouble>* distconds = ((VectorNodeObject<Superdouble>*) node.getObject(dc));
    for (int i = 0; i < nstates; i++) {
        LHOODS[i] = double(Bs.at(i) * (distconds->at(i) ));
    }
    for (int i=0; i < nstates; i++) {
        totalExp[i] = LHOODS[i];
    }
    return totalExp;
    }
}

void StateReconstructor::set_store_p_matrices(bool i) {
    store_p_matrices = i;
}

void StateReconstructor::set_use_stored_matrices(bool i) {
    use_stored_matrices = i;
}

StateReconstructor::~StateReconstructor() {

}

