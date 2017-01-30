
//#include <boost/bimap.hpp>
#include <iostream>
#include <curses.h>
#include <stdlib.h>
#include <math.h> 
using namespace std;

#include <vector>
#include <bitset>
#include <map>



typedef std::bitset<3> Sdna;
typedef unsigned int  uint;


/*
struct Comparer {
    bool operator() (const bitset<3> &b1, const bitset<3> &b2) const {
        return b1.to_ulong() < b2.to_ulong();
    }
};*/

struct base{
    Sdna code;           // store AGCTN
    uint quality = 0;    // store quality as int converted
};

//typedef vector< Sdna, char > dnamap;
//typedef vector< char, Sdna > bitmap;

typedef std::vector < base > Tdna; 


//typedef boost::bimap< boost::bimaps::set_of<Sdna *, dnasetCMP> , int> dnamap;
//typedef dnamap::value_type dnamapVT;
typedef char baseT;

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
	    void printColoredDNA(WINDOW*);




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
		
		color IntToColor(int i);
};
