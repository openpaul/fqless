
#include <iostream>
#include <curses.h>
#include <list>
//#include <string.h>
#include <stdlib.h>
#include <math.h> 
using namespace std;

#include <vector>
#include <bitset>
#include <map>
#include <future>
#include "options.h"


typedef std::bitset<3> Sdna;
typedef unsigned int  uint;



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
		
		// manipulating the sequence
		void append(string&);
		void append(Tdna&);
		void addQuality(string&);

		// retriving the sequence
		string getSequence();
	    void printColoredDNA(WINDOW*, std::pair<uint, uint> );




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
