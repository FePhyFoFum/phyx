#ifndef _BD_SIM_H_
#define _BD_SIM_H_

#include <vector>
#include <random>

class Tree; // forward declaration
class Node; // forward declaration

class BirthDeathSimulator {
private:
    int numfailures_;
    int maxfailures_;
    double birthrate_;
    double deathrate_;
    double sumrate_;
    double relative_birth_rate_;
    int extantstop_;
    double timestop_;
    int numbirth_;
    int numdeath_;
    long int seed_;
    double currenttime_;
    std::vector<Node*> extantnodes_;
    std::vector<Node*> dead_nodes_;
    std::map<Node*, double> BIRTHTIME_;
    std::map<Node*, double> DEATHTIME_;
    Node* root_;
    Tree* tree_;
    
    std::mt19937 generator_;
    std::uniform_real_distribution<double> uniformDistrib_;

    bool check_stop_conditions ();
    double time_to_next_event ();
    void event ();
    void node_death (Node *);
    void node_birth (Node *);
    void delete_dead_nodes ();
    void setup_parameters ();
    bool event_is_birth ();
    void delete_a_node (Node *);
    
    // these are not used
    double get_distance_from_tip (Node *innode);
    void set_distance_to_tip ();

public:
    BirthDeathSimulator ();
    BirthDeathSimulator (const double& estop, const double& tstop, const double& brate,
        const double& drate, const long int& seed);
    Tree * make_tree (const bool& show_dead);
    std::string get_sim_summary ();
    //~BirthDeathSimulator();
};

#endif /* _BD_SIM_H_ */
