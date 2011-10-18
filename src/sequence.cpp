/*
 * This is an extremely stripped down sequence class that is just meant to be
 * transparent and lightweight. As functionality increases, so will the 
 * complexity of the class.
 *
 */

#include <string>

#include <iostream>


using namespace std;

#include "sequence.h"

Sequence::Sequence():id(),seq(),aligned(){}

Sequence::Sequence(string _id, string _seq, bool _aligned){
    id = _id;
    seq = _seq;
    aligned = _aligned;
}
Sequence::Sequence(string _id, string _seq){
    id = _id;
    seq = _seq;
    aligned = false;
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

//this one is private
//it should eventually be generalized for both nucleotide and protein, 
//but for now it is just nucleotide for simplicity
//TODO: need to add ambiguity
string Sequence::reverse(string charin){
    string ret;
    if (charin == "-")
	ret = "-";
    if (charin == " ")
	ret = " ";
    if (charin == "A" || charin == "a"){
	ret = "T";
    }else if(charin == "T" || charin == "t"){
	ret = "A";
    }else if(charin == "C" || charin == "c"){
	ret = "G";
    }else if (charin == "G" || charin == "g"){
	ret = "C";
    }else if(charin == "U" || charin == "u"){
	ret = "A";
    }else if(charin == "m" || charin == "M"){
	ret = "K";
    }else if(charin == "r" || charin == "R"){
	ret = "Y";
    }else if(charin == "y" || charin == "Y"){
	ret = "R";
    }else if(charin == "k" || charin == "K"){
	ret = "M";
    }else if(charin == "v" || charin == "V" ){
	ret = "B";
    }else if(charin == "h" || charin == "H" ){
	ret = "D";
    }else if(charin == "d" || charin == "D" ){
	ret = "H";
    }else if(charin == "b" || charin == "B" ){
	ret = "V";
    }else if (charin == "n" || charin == "N" || charin == "x" || charin == "X"){
	ret = "N";
    }
    return ret;
}

string Sequence::reverse_complement(){
    string rcomp = seq;
    for (unsigned int i=0 ;i < rcomp.size(); i++){
	rcomp.replace(i,1,reverse(seq.substr(seq.size()-i-1,1)));
    }
    return rcomp;
}

void Sequence::perm_reverse_complement(){
    string rcomp = seq;
    for (unsigned int i=0 ;i < rcomp.size(); i++){
	rcomp.replace(i,1,reverse(seq.substr(seq.size()-i-1,1)));
    }
    seq = rcomp;
}
