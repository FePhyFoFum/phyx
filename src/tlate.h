/*
 * tlate.h
 *
 *  Created on: Jun 16, 2015
 *      Author: joe
 */

#ifndef _TLATE_H_
#define _TLATE_H_

using namespace std;

class TLATE {
public:
    TLATE();
    string Codon_to_AA (string& sequences);
    string Translate (string& dna);
    //virtual ~tlate();
};

#endif /* _TLATE_H_ */
