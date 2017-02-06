#include "fastq.h"




fastq::fastq(const char* inPath){

    string line; 
    isFastq     = false;
    isFasta     = false;
    file        = inPath;

    // determine the amount of seq in the file
    //countSequences(inPath);

    return;
}


void fastq::load2show(options* opts){
    string line;
    uint relPos = 0;
    uint i = 0;
    int number;
    string name;
    string dnaseq;
    string qualseq;
    
    // clear what we have
    content.clear();
    
    // open file
    ifstream infile;
	infile.open(file);
	
	if(opts->lastInPad == 0){
	    return;
	}
	
	// jump to start of what we want to load
	infile.seekg(index[opts->firstInPad].tellg);
	i = opts->firstInPad;
	
	 while( i < opts->lastInPad and getline(infile, line) ){

	    
	        // each sequence of fastq
	        // fills at least two lines, name and spacer
	        // then the chars are counted and divisted by cols
            
            
	        if(relPos == 0){
			    // this is the name, line
			    // make new entry into content
			    name     = line.erase(0, 1);;
			    number   = i;
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
	            b.inpad         = true;
	            setDNAline(b, dnaseq);
	            addQualityData(b,qualseq);
	            

	            // now save it
                content.push_back(b);
	           
	           
		    }
		    
			relPos++;
			if( relPos == 4 ){
			    // reset the relative Position counter
			    i++;
	            relPos = 0;
			}
			

	    }
	
	
}

int fastq::readmore(options* opts, int dir, WINDOW* Wtext){
    //dir 1 down 0 up
    string line;
    uint bfull;
    bfull       = 0;
    uint i      = content.size();
    uint relPos = 0;
    
    
    ifstream infile;
	infile.open(file);
	
	
	if(infile.is_open() ){
	
	    if( dir == 1  ){
            // jump down
            infile.seekg(opts->tellg);
        }
	
	    int number;
        string name;
        string dnaseq;
        string qualseq;
        double ctellg;
            
	    while( bfull < opts->avaiLines and getline(infile, line) ){

	    
	        // each sequence of fastq
	        // fills at least two lines, name and spacer
	        // then the chars are counted and divisted by cols
            
            
	        if(relPos == 0){
			    // this is the name, line
			    // make new entry into content
			    name     = line.erase(0, 1);;
			    number   = i;
			    ctellg   = infile.tellg();          // remember the position of the name start
			    ctellg   = ctellg - line.length();   // tahts why we go back the name line length
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
	            b.tellg         = ctellg;
	            b.inpad         = false;
	            setDNAline(b, dnaseq);
	            addQualityData(b,qualseq);
	            
	            // update the bfull variable
	            bfull = bfull + 2 + ceil((float)b.dna.sequence.size()/(float)opts->textcols) ;
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
			
		    // save seek pointer for later readings
	        opts->tellg = infile.tellg();
	    }

	    
	    
	}

	// close up
	infile.close();
    
    return bfull;
}


void fastq::showthese(options* opts, int dir, WINDOW* Wtext){
     //dir 1 means down, 0 means up

    // sets the in pad variable to indictae which entries to show
    // also determines the offset the editor must have to make scrolling look 
    // nice
    
    // first iterate through and make a vector of lines/per entry
    // remember also which sequence was the last the user had in the pad
    // this is directionally dependent, so either the first row of the last
    // pad or the last row.
    
    vector<int> linevec; // hold the number of lines, each entry has
    int lines;
    uint lineSum    = 0;
    int firstSeq    = opts->firstInPad;
    int lastSeq     = opts->lastInPad;
    int i           = 0;
    
   
    // build a vector holding line counts
    for(auto& it: index) {
        indexStruc& ind = it;
        lines = 2 + ceil((float)ind.lengthName/(float)opts->textcols) ;
        linevec.push_back(lines);
     //   if(ind.inpad and firstSeq == -1){
       //     firstSeq = i;
       // }
       // if(ind.inpad){
       //     lastSeq = i;
       // }
        // reset inpad values, we set them later again
        //ind.inpad = false;
        //i++;
    }
    
    //if(firstSeq == -1){firstSeq = 0;}
    //if(lastSeq == -1){lastSeq = -0;}
    //wprintw(Wtext, "first %i and last %i \n",  firstSeq, lastSeq);
    //if(firstSeq == -1){firstSeq = -0;}
    
    // ok, indexing is done, now lets find out what to show
    int linesBelowAndUp = round(opts->avaiLines/2);
    if(linesBelowAndUp < opts->textrows){
        // this fixes a bug, whne the window is larger then half
        // the buffer which then allows escaping of the window
        linesBelowAndUp = opts->textrows + 1;
    }

    
    // scroll N lines up in both cases
    uint j = 0;
    if(dir == 1){
        i = lastSeq;
    }else{
        i = firstSeq;
    }
    while(j < linesBelowAndUp && i > 0){
        j = j + linevec[i];
        i--;
    }
    //wprintw(Wtext, "i %i and j %i \n",  i, j);
    // set offset to match the new pad value
    opts->offset = j;
    if(dir == 1){
        opts->offset -= opts->textrows;
    }
    if(opts->offset < 0) opts->offset = 0;

    // now we move down find the last entry in the pad
    
    //vector<uint> pos2load;
    if(i < 0){ i = 0;}
    opts->firstInPad = i;
    
    j = 0;
    while(j < opts->avaiLines && i < index.size()){
        //pos2load.push_back(i);
        j = j + linevec[i];
        i++;
        opts->lastInPad = i;
    }
    

    opts->linesTohave = j;
    wprintw(Wtext, "first %i and last %i  offset %i seq in %i\n",  opts->firstInPad, opts->lastInPad, opts->offset,  opts->lastInPad - opts->firstInPad);
    // now we know what to show, we can actually go ahead and load them.
    load2show(opts);
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
	         if( currLine % 4 == 0 and line.length()>0){
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
    infile.close();
    return;
}



void fastq::buildIndex(options* opts){
    //dir 1 down 0 up
    string line;
    uint bfull  = 0;
    uint i      = index.size();
    uint relPos = 0;
    
    // check if this is the first time we do this
    uint firstIndex = false;
    if(index.size() == 0){
        firstIndex = true;
    }
    
    
    ifstream infile;
	infile.open(file);
	
	
	if(infile.is_open() ){
	
	    if( firstIndex == false  ){
            // jump down
            // last index pos:
            infile.seekg(opts->IndexTellg);
        }
	
        uint number;
        double ctellg;
        uint lengthName;
        uint lengthSeq;
        
        
            
	    while( getline(infile, line) ){
            
            if(firstIndex and bfull >= 2*opts->avaiLines ){
                // break if we only want to index a first quick glance of the file
                // for fast response times
                break;
            }
            
	    
	        // each sequence of fastq
	        // fills at least two lines, name and spacer
	        // then the chars are counted and divisted by cols
            
            
	        if(relPos == 0){
			    // this is the name, line
			    // make new entry into content
			    lengthName      = line.erase(0, 1).size();
			    number          = i;
			    ctellg          = infile.tellg();          // remember the position of the name start
			    ctellg          = ctellg - line.size() - 2;  // tahts why we go back the name line length
			    i++;
		    }else if(relPos == 1){
			    // this should be DNA
			    lengthSeq = line.size();
			    
		    }else if(relPos == 3){
		       		    
			    // now create the fastq seq and put it in the right place
	            indexStruc ind;
	            	            
	            ind.number          = number;
	            ind.lengthSeq       = lengthSeq;
	            ind.lengthName      = lengthName;
	            ind.tellg           = ctellg;     
	            ind.inpad           = false;      
	            // update the bfull variable
	            bfull = bfull + 2 + ceil((float)lengthSeq/(float)opts->textcols) ;
	            //            name   rows of DNA              space
	            index.push_back(ind);

		    }
		    
			relPos++;
			if( relPos == 4 ){
			     // reset the relative Position counter
	            relPos          = 0;
			}
			
		    // save seek pointer for later readings
	        opts->IndexTellg = infile.tellg();
	    }   
	    
	}

	// close up
	infile.close();
    
    return;
}



void fastq::setDNAline(fastqSeq& b, string& sDNA){
	b.dna.append(sDNA);
}

void fastq::addQualityData(fastqSeq& b, string& qual){
	b.dna.addQuality(qual);   
}
