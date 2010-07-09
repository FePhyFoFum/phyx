
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;

#include "bd_sim.h"
#include "tree.h"
#include "node.h"

BirthDeathSimulator::BirthDeathSimulator(double estop,double tstop, double brate, double drate):failures(0),
		maxfailures(1000),birthrate(brate),deathrate(drate),sumrate(brate+drate),relative_birth_rate(brate/(brate+drate)),
		extantstop(estop),timestop(tstop),numofchanges(0),currenttime(0.0),extantnodes(vector<Node*>()),
		BIRTHTIME(map<Node*,double>()),DEATHTIME(map<Node*,double>()){}

BirthDeathSimulator::BirthDeathSimulator():failures(0),maxfailures(1000),birthrate(0.1),deathrate(0.05),
		sumrate(0.1+0.05),relative_birth_rate(0.1/(0.1+0.05)),extantstop(10),timestop(0),numofchanges(0),
		currenttime(0.0),extantnodes(vector<Node*>()),BIRTHTIME(map<Node*,double>()),DEATHTIME(map<Node*,double>()){}

void BirthDeathSimulator::setup_parameters(){
	numofchanges = 0;
	currenttime = 0.0;
	extantnodes = vector<Node*>();
	BIRTHTIME = map<Node*,double>();
	DEATHTIME = map<Node*,double>();
}

Tree * BirthDeathSimulator::make_tree(bool show_dead){
	setup_parameters();
	root = new Node();
	BIRTHTIME[root] = currenttime;
	extantnodes.push_back(root);
	while(check_stop_conditions()){
		double dt = time_to_next_sp_event();
		currenttime += dt;
		event();
		if(extantnodes.size() < 1){
			failures += 1;
			setup_parameters();
			root = new Node();//need to clean
			BIRTHTIME[root] = currenttime;//need to clean
			extantnodes.push_back(root);//empty
		}
	}
	vector<Node*> temp_extant_nodes(extantnodes);
	for (unsigned int i=0;i<temp_extant_nodes.size();i++){
		try{
			DEATHTIME[temp_extant_nodes[i]];
		}catch( char * str ){
			node_death(temp_extant_nodes[i]);
			//temp_extant_nodes[i]->istip = 1
		}
	}
	if(show_dead == false){
		delete_dead_nodes();
	}
	root->setBL(0);
	double totallength = 0;
	int count = 1;
	Tree * tree = new Tree(root);
	for(int i=0;i<tree->getExternalNodeCount();i++){
		totallength += tree->getExternalNode(i)->getBL();
		std::stringstream out;
		out << count;
		tree->getExternalNode(i)->setName("taxon_"+out.str());
		count += 1;
	}
	return tree;
}

bool BirthDeathSimulator::check_stop_conditions(){
	bool keepgoing = true;
	if (extantstop > 0){
		if (extantnodes.size() >= extantstop){
			currenttime += time_to_next_sp_event();
			keepgoing = false;
		}
	}
	if(timestop > 0){
		if (currenttime >= timestop){
			currenttime = timestop;
			keepgoing = false;
		}
	}
	return keepgoing;
}

double BirthDeathSimulator::time_to_next_sp_event(){
	return  0.0;
}

void BirthDeathSimulator::event(){

}

void BirthDeathSimulator::node_death(Node *innode){

}

void BirthDeathSimulator::delete_dead_nodes(){

}
