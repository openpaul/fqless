#include <iostream>
#include <fstream>
//#include <vector>
//#include <bitset>
#include <getopt.h> 
#include <string> 
#include <stdlib.h> 
using namespace std;

#include "DNA.h"
const int SEQLIMIT = 500; // currently i limit this here, to prevent crashs, until a better mechanism is implemented


struct fastqSeq{
    DNA dna;
    string name;
    uint number; // is this the first, seocnd or third sequence of the file.
};

class fastq{
    
    public:


        //string name;    
        vector<fastqSeq> content;

        
        int nOfSequences = 0;
        
        // read fastq file from disk
        fastq(const char*);   
    
        

	private:
	    uint currentSequence = 0; // keep track of which sequence is currently open, read, printed
		void setDNAline(string&); 
		void addQualityData(string&);
		
		string quality = "sanger"; // encoding of the file quality string etc.
		                           // may be  sanger, solexa, illumina1.3, illumina1.5 illumina1.8 
		                           // for color chemes required
        void countSequences(const char*);
        
        bool isFastq;
        bool isFasta;

};
