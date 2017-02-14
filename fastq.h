
#ifndef FASTQ_H
#define FASTQ_H

#include <iostream>
#include <fstream>
//#include <vector>
//#include <bitset>
#include <getopt.h> 
#include <string> 
#include <stdlib.h> 
using namespace std;

#include "DNA.h"






struct fastqSeq{
    DNA dna;
    string name;
    
    uint number; // is this the first, seocnd or third sequence of the file.
    double tellg; // index position in the source file
    bool inpad;
};


struct indexStruc{
    uint number;
    double tellg;
    uint lengthName;
    uint lengthSeq; 
    bool inpad;
};

class fastq{
    
    public:


        //string name;    
        vector<fastqSeq> content;
        vector<indexStruc> index;

        
        int nOfSequences = 0;
        
        
        int seqOnScreen;
        
        // read fastq file from disk
        fastq(string);  
        // read n terminal lines from file, manage which line was the last read one (tellg)
        fastq(const char*, int lines, int cols); 
        
        int readmore(options*, int , WINDOW*);
        void showthese(options* , int ,WINDOW*);
        
        void buildIndex(options*);
        void load2show(options*);
        
        // posible quality values
        vector<string> possibleQual;
    
        string file;

        int minQal =  100000;
        int maxQal = -100000;
	private:
	    uint currentSequence = 0; // keep track of which sequence is currently open, read, printed
		void setDNAline(fastqSeq& b, string&); 
		void addQualityData(fastqSeq& b, string&, options*);
		
		string quality = "sanger"; // encoding of the file quality string etc.
		                           // may be  sanger, solexa, illumina1.3, illumina1.5 illumina1.8 
		                           // for color chemes required
        void countSequences(const char*);
        



        
        

        
        
        

};
#endif
