#ifndef DNA_H
#define DNA_H

#include "options.h"

#include <iostream>
#include <stdio.h>
#include <list>
#include <vector>
#include <bitset>
#include <map>



#include <curses.h>
using namespace std;

typedef char Sdna;



struct base{
    Sdna code;           // store AGCTN
    uint quality = 0;    // store quality as int converted
};



typedef std::vector < base > Tdna; 

struct color {
    uint R;
    uint G;
    uint B;
};




class DNA{

public:

    Tdna getSeq(){return(sequence);}
    Tdna sequence; // we use three bits to store AGCT and possible N X R etc..

    // create new DNA element
    DNA(vector < char >);
    DNA();

    uint size(){
        return(sequence.size());
    };

    // manipulating the sequence
    void append(string&);
    void append(Tdna&);
    void addQuality(string&);

    // retriving the sequence
    string getSequence();
    void printColoredDNA(WINDOW*, std::pair<uint, uint>, bool );




private:
    // convert Character to bit and reverse
    // this should save memory, as one Char only takes 3 instead of 8 bits to store.
    base char2DNA(char &);
    char DNA2char(base &);

    // DNA and bit mapping values and functions
    //map <char, Sdna> getDNAMap(); 	
    //bitmap getBitMap();
    //dnamap getDNAMap();
    //dnamap dnaMap = getDNAMap() ;

    //		map <char, Sdna> dnaMap = getDNAMap();	
    //bitmap bitMap = getBitMap();

    color IntToColor(int i, std::pair<uint, uint> );
};
#endif
