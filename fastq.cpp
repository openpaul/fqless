
#include "fastq.h"
#include <math.h> 

fastq::fastq(string inPath){

    file        = inPath;

    return;
}


void fastq::load2show(options* opts){
    string line;
    uint relPos    = 0;
    uint i         = 0;
    int number     = 0;
    string name    = "";
    string dnaseq  = "";
    string qualseq = "";

    // clear what we have
    content.clear();

    // open file
    gzFile infile = gzopen(opts->input, "rb");
    //ifstream infile;
    //infile.open(opts->input);

    if(opts->lastInPad == 0){
        return;
    }

    // jump to start of what we want to load
    gzseek(infile, index[opts->firstInPad].tellg, SEEK_SET);
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
            addQualityData(b,qualseq, opts);


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
    int firstSeq    = opts->firstInPad;
    int lastSeq     = opts->lastInPad;
    int i           = 0;


    // build a vector holding line counts
    for(auto& it: index) {
        indexStruc& ind = it;
        lines = 1 + ceil((float)ind.lengthName/(float)opts->textcols) +
        ceil((float)ind.lengthSeq/(float)opts->textcols) ;


        linevec.push_back(lines);
    }



    // ok, indexing is done, now lets find out what to show
    int linesBelowAndUp = round(opts->avaiLines/2);
    if(linesBelowAndUp < opts->textrows){
        // this fixes a bug, whne the window is larger then half
        // the buffer which then allows escaping of the window
        linesBelowAndUp = opts->textrows + 1;
    }


    // scroll N lines up in both cases
    int j = 0;
    if(dir == 1){
        i = lastSeq;
    }else{
        i = firstSeq;
    }

    while(j < linesBelowAndUp && i > 0){
        j = j + linevec[i];
        i--;
    }


    // set offset to match the new pad value
    opts->offset = j;
    if(dir == 1){
        opts->offset -= opts->textrows;
    }
    if(opts->offset < 0) opts->offset = 0;

    // now we move down find the last entry in the pad
    if(i < 0){ i = 0;}
    opts->firstInPad = i;

    j = 0;
    while(j < (int)opts->avaiLines && i < (int)index.size()){
        j = j + linevec[i];
        i++;
        opts->lastInPad = i;
    }

    opts->linesTohave = j;
    // now we know what to show, we can actually go ahead and load them.
    load2show(opts);
}


int fastq::getline(gzFile& infile, string& line){
    line = "";
    int a;
    char aa;
    while(true){
        a =  gzgetc(infile);
        if(a == -1){
            return(0);
        }
        aa = (char)a;
        if(aa == '\n' or aa == '\r'){
            break;
        }
        line.push_back(aa);
    }
    return 1;
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


    gzFile infile = gzopen(opts->input, "rb");


    if(infile != NULL ){

        if( firstIndex == false  ){
            // jump down
            // last index pos:
            gzseek(infile, opts->IndexTellg, SEEK_SET);
        }

        uint number         = 0;
        double ctellg       = 0;
        uint lengthName     = 0;
        uint lengthSeq      = 0;



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
                ctellg          = gztell(infile);          // remember the position of the name start
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
            opts->IndexTellg = gztell(infile);
        }   

    }

    // close up
    gzclose(infile);

    return;
}



void fastq::setDNAline(fastqSeq& b, string& sDNA){
    b.dna.append(sDNA);
}

void fastq::addQualityData(fastqSeq& b, string& qual, options* opts){
    b.dna.addQuality(qual);

    //estimate quality range
    int q;
    for(auto it = qual.begin(); it != qual.end(); ++it) {
        q = static_cast<int>(*it);
        if(q < minQal){
            minQal = q;
        }
        if(q > maxQal){
            maxQal = q;
        }
    }
    // set possible quality
    possibleQual.clear();
    for(auto& qv : opts->qm) {
        if(minQal >= qv.second.first && maxQal <= qv.second.second){
            possibleQual.push_back(qv.first);
        }
    }
    // in case there is no kown mode
    // the integers still code for quality
    // so make a new unknown coding
    if(possibleQual.size() == 0){
        opts->qm["unknown"] = std::make_pair(minQal, maxQal);
        possibleQual.push_back("unknown");
    }

}

qualitymap buildQualityMap(){
    qualitymap map;
    //   Name                     //min  max
    map["Sanger"]           = std::make_pair(33, 73);
    map["Solexa"]           = std::make_pair(59, 104);
    map["Illumina 1.3+"]    = std::make_pair(64, 104);
    map["Illumina 1.8+"]    = std::make_pair(33, 74);
    return map;
}
