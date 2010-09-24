

#include <string>

using namespace std;

#include "sequence.h"

Sequence::Sequence():id(),seq(),aligned(){}

Sequence::Sequence(string _id, string _seq, bool _aligned){
	id = _id;
	seq = _seq;
	aligned = _aligned;
}

bool Sequence::is_aligned(){
	return aligned;
}

string Sequence::get_sequence(){
	return seq;
}

string Sequence::get_id(){
	return id;
}

void Sequence::set_sequence(string _seq){
	seq = _seq;
}

void Sequence::set_id(string _id){
	id = _id;
}

void Sequence::set_aligned(bool _aligned){
	aligned = _aligned;
}
