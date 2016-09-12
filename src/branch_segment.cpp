/*
 * BranchSegment.cpp
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 */


#include "branch_segment.h"
#include "rate_model.h"

#include <vector>
using namespace std;

BranchSegment::BranchSegment(double dur,int per):duration(dur),period(per),
        model(NULL),fossilareaindices(vector<int>()),startdistint(-666),
        distconds(NULL),ancdistconds(NULL) {}

void BranchSegment::setModel(RateModel * mod) {
    model = mod;
}

/*void BranchSegment::setStartDist(vector<int> sd) {
    startdist = sd;
}*/

void BranchSegment::clearStartDist() {
    //startdist.clear();
    startdistint = -666; //null is -666
}

double BranchSegment::getDuration() {
    return duration;
}

int BranchSegment::getPeriod() {
    return period;
}
/*
vector<int> BranchSegment::getStartDist() {
    return startdist;
}*/

void BranchSegment::set_start_dist_int(int d) {
    startdistint = d;
}

int BranchSegment::get_start_dist_int() {
    return startdistint;
}

RateModel * BranchSegment::getModel() {
    return model;
}

vector<int> BranchSegment::getFossilAreas() {
    return fossilareaindices;
}

void BranchSegment::setFossilArea(int area) {
    fossilareaindices.push_back(area);
}
