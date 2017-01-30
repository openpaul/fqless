#include "fastq.h"

fastq::fastq(const char*  inPath){
    
	string line; 
    isFastq = false;
    isFasta = false;
    
    // determine the amount of seq in the file
    countSequences(inPath);
    
	if(inPath == NULL){
		return;
	}
	ifstream infile;
	infile.open(inPath);

	if(infile.is_open()){
		
		uint i = 0;
		uint relPos = 0;

		while( getline(infile, line) ){
			if( isFastq and i > 4 * SEQLIMIT) break;
			if(relPos == 0){
				// this is the name, line
				// make new entry into content
				fastqSeq b;
				content.push_back(b);
				content.back().name     = line;
				content.back().number   = i;

			}
			
			if(relPos == 1){
				// this should be DNA
				setDNAline(line);
			}
			if(relPos == 2){
				
				//setDNAline(line);
			}
				
			// keep track of relative position
			// this way we always know what kind of line
			// to expect
			if(relPos == 3){
			    // this should quality data
                addQualityData(line);
				relPos = 0;
			}else{
				relPos++;
			}
			i++; // count sequences
		}
		
		infile.close();
		
		//cout << "read " << content.size() << " sequences" << std::endl;
	}else{
		cout << "File is not existing or unable to open" << std::endl;
	}
}

void fastq::countSequences(const char* inPath){
    /*
        this function actulaly opens the file and has a look at it
        it can decide if the file is fasta or fastq
        also it counts the number of sequences stored in the file
    */
    string line; 

	if(inPath == NULL){
		return;
	}
	ifstream infile;
	infile.open(inPath);
        
    int currLine = 0;
	if(infile.is_open()){
	    while( getline(infile, line) ){
	        if( currLine % 4 == 0 ){
	            char d;
	            d = line.at(0);
	             if( currLine == 0){
	                if( d == '@'){ isFastq = TRUE; }
	                if( d == '>'){ isFasta = TRUE; }
	             }
	            if( ( isFastq and d == '@' ) or( isFasta and d == '>' ) ){
	                nOfSequences++;
	            }
	            
	        }
	        currLine++ ;
	    }
    }
    return;
}


void fastq::setDNAline(string& sDNA){
	content.back().dna.append(sDNA);
}

void fastq::addQualityData(string& qual){
	content.back().dna.addQuality(qual);   
}
