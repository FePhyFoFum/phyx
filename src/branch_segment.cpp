#include <vector>

#include "branch_segment.h"
#include "rate_model.h"


BranchSegment::BranchSegment(double dur, int per):duration(dur), period(per),
        model(nullptr), fossilareaindices(std::vector<int>()), startdistint(-666),
        distconds(nullptr), ancdistconds(nullptr) {}


void BranchSegment::setModel (RateModel * mod) {
    model = mod;
}


/*void BranchSegment::setStartDist (std::vector<int> sd) {
    startdist = sd;
}*/


void BranchSegment::clearStartDist () {
    //startdist.clear();
    set_start_dist_int(-666); //null is -666
}


double BranchSegment::getDuration () const {
    return duration;
}


int BranchSegment::getPeriod () const {
    return period;
}


/*
vector<int> BranchSegment::getStartDist () const {
    return startdist;
}*/


void BranchSegment::set_start_dist_int (int d) {
    startdistint = d;
}


int BranchSegment::get_start_dist_int () const {
    return startdistint;
}


RateModel * BranchSegment::getModel () const {
    return model;
}


std::vector<int> BranchSegment::getFossilAreas () const {
    return fossilareaindices;
}


// not used
void BranchSegment::setFossilArea (int area) {
    fossilareaindices.push_back(area);
}
