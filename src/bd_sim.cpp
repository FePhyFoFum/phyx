#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <random>

#include "bd_sim.h"
#include "tree.h"
#include "node.h"
#include "utils.h"


BirthDeathSimulator::BirthDeathSimulator (const double& estop, const double& tstop,
    const double& brate, const double& drate, const int& seed):failures_(0), maxfailures_(1000), birthrate_(brate),
    deathrate_(drate), sumrate_(brate+drate), relative_birth_rate_(brate/(brate+drate)),
    extantstop_(estop), timestop_(tstop), numofchanges_(0), currenttime_(0.0),
    extantnodes_(std::vector<Node*>()), BIRTHTIME_(std::map<Node*, double>()),
    DEATHTIME_(std::map<Node*, double>()) {
    if (seed == -1) {
        generator_ =  std::mt19937(get_clock_seed());
    } else {
        generator_ =  std::mt19937(seed);
    }
    uniformDistrib_ =  std::uniform_real_distribution<double>(0.0, 1.0);
}


// not currently used
BirthDeathSimulator::BirthDeathSimulator ():failures_(0), maxfailures_(1000),
    birthrate_(0.1), deathrate_(0.05), sumrate_(0.1+0.05),
    relative_birth_rate_(0.1/(0.1+0.05)), extantstop_(10), timestop_(0), numofchanges_(0),
    currenttime_(0.0), extantnodes_(std::vector<Node*>()), BIRTHTIME_(std::map<Node*, double>()),
    DEATHTIME_(std::map<Node*, double>()) {
        generator_ =  std::mt19937(get_clock_seed());
        uniformDistrib_ =  std::uniform_real_distribution<double>(0.0, 1.0);
}


void BirthDeathSimulator::setup_parameters () {
    numofchanges_ = 0;
    currenttime_ = 0.0;
    extantnodes_ = std::vector<Node*>();
    dead_nodes_ = std::vector<Node*>();
    BIRTHTIME_ = std::map<Node*, double>();
    DEATHTIME_ = std::map<Node*, double>();
}


Tree * BirthDeathSimulator::make_tree (const bool& show_dead) {
    setup_parameters();
    root_ = new Node();
    BIRTHTIME_[root_] = currenttime_;
    extantnodes_.push_back(root_);
    
    // actually want to start with 2 lineages
    node_birth(extantnodes_[0]);
    
    // reset failures to zero. don't want to accumulate errors across replicates
    failures_ = 0;
    bool going = true;
    while (going) {
        double dt = time_to_next_event();
        currenttime_ += dt;
        going = check_stop_conditions();
        if (going) {
            event();
            if (extantnodes_.size() < 1) {
                failures_ += 1;
                //std::cout << "failed!" << std::endl;
                if (failures_ >= maxfailures_) {
                    std::cout << "Reached maximum number of failures (" << failures_
                        << "). Quitting." << std::endl;
                    exit(0);
                }
                setup_parameters();
                root_ = new Node();//need to clean
                BIRTHTIME_[root_] = currenttime_;//need to clean
                extantnodes_.push_back(root_);//empty
            }
        }
    }
    
    std::vector<Node*> temp_extant_nodes(extantnodes_);
    for (unsigned int i = 0; i < temp_extant_nodes.size(); i++) {
        try {
            DEATHTIME_[temp_extant_nodes[i]];
            node_death(temp_extant_nodes[i]);
        } catch( char * str ) {
            std::cout << "catch" << std::endl;
            node_death(temp_extant_nodes[i]);
            //temp_extant_nodes[i]->istip = 1
        }
    }

    root_->setBL(0);
    double totallength = 0;
    int count = 1;
    tree_ = new Tree(root_);
    tree_->processRoot();
    if (show_dead == false) {
        delete_dead_nodes();
    }

    for (int i = 0; i < tree_->getExternalNodeCount(); i++) {
        totallength += tree_->getExternalNode(i)->getBL();
        std::stringstream out;
        out << count;
        tree_->getExternalNode(i)->setName("taxon_"+out.str());
        count += 1;
    }
    return tree_;
}


bool BirthDeathSimulator::check_stop_conditions () {
    bool keepgoing = true;
    if (extantstop_ > 0) {
        if ((int)extantnodes_.size() >= extantstop_) {
            // this ensures tips do not have 0 edge lengths
            currenttime_ += time_to_next_event();
            keepgoing = false;
        }
    }
    if (timestop_ > 0) {
        if (currenttime_ >= timestop_) {
            // definitely want this, or tree can be older than timestop
            currenttime_ = timestop_;
            keepgoing = false;
        }
    }
    return keepgoing;
}


// time until next event (rate = sum of birth + death)
double BirthDeathSimulator::time_to_next_event () {
    //double num = rand() / double(RAND_MAX);
    double num = uniformDistrib_(generator_);
    return (-log(num))/ ( (extantnodes_.size()) * sumrate_);
}


void BirthDeathSimulator::event () {
    std::uniform_int_distribution<int> intDistrib(0, (extantnodes_.size()-1));
    int random_integer = intDistrib(generator_);
    //std::cout << "extantnodes.size() = " << extantnodes.size() << "; random_integer = "
    //    << random_integer << std::endl;
    Node * extant = extantnodes_[random_integer];
    if (event_is_birth()) {
        node_birth(extant); // real speciation
    } else {
        node_death(extant);
        dead_nodes_.push_back(extant);
    }
}


void BirthDeathSimulator::node_death (Node *innode) {
    DEATHTIME_[innode] = currenttime_;
    double bl = DEATHTIME_[innode] - BIRTHTIME_[innode];
    innode->setBL(bl);
    extantnodes_.erase(find(extantnodes_.begin(), extantnodes_.end(), innode));
}


void BirthDeathSimulator::node_birth (Node *innode) {
    Node * left = new Node();
    Node * right = new Node();
    innode->addChild(*left);
    innode->addChild(*right);
    BIRTHTIME_[left] = currenttime_;
    BIRTHTIME_[right] = currenttime_;
    node_death(innode);
    extantnodes_.push_back(left);
    extantnodes_.push_back(right);
}


void BirthDeathSimulator::delete_dead_nodes () {
    for (unsigned int i = 0; i < dead_nodes_.size(); i++) {
        delete_a_node(dead_nodes_[i]);
    }
}


void BirthDeathSimulator::set_distance_to_tip () {
    for (int i = 0; i < tree_->getExternalNodeCount(); i++) {
        double curh = 0.0;
        tree_->getExternalNode(i)->setHeight(curh);
        Node * tnode = tree_->getExternalNode(i);
        while (tnode->hasParent()) {
            curh += tnode->getBL();
            if (tnode->getHeight() < curh) {
                tnode->setHeight(curh);
            }
            tnode = tnode->getParent();
        }
        curh += tnode->getBL();
        if (tnode->getHeight()<curh) {
            tnode->setHeight(curh);
        }
    }
}


void BirthDeathSimulator::delete_a_node (Node * innode) {
    Node * tparent = innode->getParent();
    if (tparent != root_) {
        Node * child = NULL;
        for (int i = 0; i < tparent->getChildCount(); i++)
            if (tparent->getChild(i) != innode) {
                child = tparent->getChild(i);
            }
        Node * pparent = tparent->getParent();
        tparent->removeChild(*innode);
        tparent->removeChild(*child);
        pparent->removeChild(*tparent);
        pparent->addChild(*child);
        child->setParent(*pparent);
        child->setBL(child->getBL()+tparent->getBL());
    } else {
        Node * child = NULL;
        for (int i = 0; i < tparent->getChildCount(); i++)
            if (tparent->getChild(i) != innode) {
                child = tparent->getChild(i);
            }
        tparent->removeChild(*innode);
        tree_->setRoot(child);
        root_ = child;
    }
}


bool BirthDeathSimulator::event_is_birth () {
    double x = uniformDistrib_(generator_);
    if (x < relative_birth_rate_) {
        return true;
    } else {
        return false;
    }
}


double BirthDeathSimulator::get_distance_from_tip (Node *innode) {
    Node * cur = innode;
    double curh = 0.0;
    while (cur->hasParent()) {
        curh += cur->getBL();
        cur = cur->getParent();
    }
    curh += cur->getBL();
    return curh;
}
