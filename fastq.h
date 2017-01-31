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

const int nSCREENS = 20; // Number of screen to store additionally to the existing screen. 



struct fastqSeq{
    DNA dna;
    string name;
    uint number; // is this the first, seocnd or third sequence of the file.
};

class fastq{
    
    public:


        //string name;    
        list<fastqSeq> content;

        
        int nOfSequences = 0;
        
        
        int seqOnScreen;
        
        // read fastq file from disk
        fastq(const char*);  
        // read n terminal lines from file, manage which line was the last read one (tellg)
        fastq(const char*, int lines, int cols); 
        
        int readmore(int, int , int , int , WINDOW*);
        
    
        

	private:
	    uint currentSequence = 0; // keep track of which sequence is currently open, read, printed
		void setDNAline(fastqSeq& b, string&); 
		void addQualityData(fastqSeq& b, string&);
		
		string quality = "sanger"; // encoding of the file quality string etc.
		                           // may be  sanger, solexa, illumina1.3, illumina1.5 illumina1.8 
		                           // for color chemes required
        void countSequences(const char*);
        
        const char* file;
        bool isFastq;
        bool isFasta;
        

};
