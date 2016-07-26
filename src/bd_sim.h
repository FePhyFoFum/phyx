#ifndef _BD_SIM_H_
#define _BD_SIM_H_

#include <random>

using namespace std;

#include "tree.h"
#include "node.h"

class BirthDeathSimulator {
private:
    int failures_;
    int maxfailures_;
    double birthrate_;
    double deathrate_;
    double sumrate_;
    double relative_birth_rate_;
    int extantstop_;
    double timestop_;
    int numofchanges_;
    double currenttime_;
    vector<Node*> extantnodes_;
    vector<Node*> dead_nodes_;
    map<Node*,double> BIRTHTIME_;
    map<Node*,double> DEATHTIME_;
    Node* root_;
    Tree* tree_;
    
    mt19937 generator_;
    std::uniform_real_distribution<double> uniformDistrib_;

    bool check_stop_conditions();
    double time_to_next_event();
    void event();
    void node_death(Node *);
    void node_birth(Node *);
    void delete_dead_nodes();
    void setup_parameters();
    bool event_is_birth();
    void delete_a_node(Node *);
    double get_distance_from_tip(Node *innode);
    void set_distance_to_tip();

public:
    BirthDeathSimulator();
    BirthDeathSimulator(double estop, double tstop, double brate, double drate, int seed);

    Tree * make_tree(bool);

    //~BirthDeathSimulator();
};

#endif /* _BD_SIM_H_ */
