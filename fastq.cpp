#include "fastq.h"




fastq::fastq(const char* inPath){

    string line; 
    isFastq     = false;
    isFasta     = false;
    file        = inPath;

    // determine the amount of seq in the file
    countSequences(inPath);

    return;
}

int fastq::readmore(int offset, int bsize, int cols, int dir, WINDOW* Wtext){
    //dir 1 down 0 up
    string line;
    uint bfull;
    bfull = 0;
    uint i = 0;
    uint relPos = 0;
    

    
    if( offset != 0 ){
        // need to skip some rows apparently
    }
    
    ifstream infile;
	infile.open(file);
	if(infile.is_open()){
	
	    int number;
        string name;
        string dnaseq;
        string qualseq;
            
	    while( bfull < bsize and getline(infile, line)){
	    
	        // each sequence of fastq
	        // fills at least two lines, name and spacer
	        // then the chars are counted and divisted by cols
            
            
	        if(relPos == 0){
			    // this is the name, line
			    // make new entry into content
			    name     = line.erase(0, 1);;
			    number   = i;
			    i++;
		    }else if(relPos == 1){
			    // this should be DNA
			    dnaseq = line;
			    
		    }else if(relPos == 3){
		       
		      
		        // this should quality data
                qualseq         = line;
			    
			    // now create the fastq seq and put it in the right place
	            fastqSeq b;
	            	            
	            b.name          = name;
	            b.number        = number;
	            setDNAline(b, dnaseq);
	            addQualityData(b,qualseq);
	            
	            // update the bfull variable
	            bfull = bfull + 2 + ceil((float)b.dna.sequence.size()/(float)cols) ;
	            //            name   rows of DNA              space

	            // now put the sequence into the content variable, so it can be accessed
	            if(dir == 1){
	                content.push_back(b);
	            }else{
	                // else put each sequence at a certain position 
	                // in front of the content we already have
	                //content.push_back(b);
	                auto it = content.begin();
	                advance(it, i);
	                content.insert(it, b);
	            }
	           
		    }
		    
			relPos++;
			if( relPos == 4 ){
			     // reset the relative Position counter
	            relPos          = 0;
			}
		    
	    }
	}
	
	// close up
	infile.close();
    
    return bfull;
}

/*
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
}*/

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


void fastq::setDNAline(fastqSeq& b, string& sDNA){
	b.dna.append(sDNA);
}

void fastq::addQualityData(fastqSeq& b, string& qual){
	b.dna.addQuality(qual);   
}
