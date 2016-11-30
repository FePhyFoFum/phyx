/*
 * tlate.h
 *
 *  Created on: Jun 16, 2015
 *      Author: joe
 */

#ifndef _TLATE_H_
#define _TLATE_H_

#include <map>

using namespace std;

class TLATE {
private:
    map <string, string> table_;
    
public:
    TLATE (string const& table);
    string translate (string& dna);
    //virtual ~tlate();
};

#endif /* _TLATE_H_ */
