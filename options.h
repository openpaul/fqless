
#ifndef OPTIONS_H
#define OPTIONS_H
typedef std::map<std::string, std::pair<int, int> > qualitymap;


static qualitymap buildQualityMap(){
    qualitymap map;
    //   Name                     //min  max
    map["Sanger"]           = std::make_pair(33, 73);
    map["Solexa"]           = std::make_pair(59, 104);
    map["Illumina 1.3+"]    = std::make_pair(64, 104);
    map["Illumina 1.8+"]    = std::make_pair(33, 74);
    return map;
}




struct options{
    string input;
    
    bool linenumbers;
    int buffersize; 
   
    
     
    
    // actually not options but the current state of the software
    int textrows;
    int textcols;
    int offset;
    
    int qualitycode;
   
    uint avaiLines;
    uint linesTohave;
    
 // reading the file
    unsigned long int tellg; 
    unsigned long int IndexTellg;
    
    uint firstInPad;
    uint lastInPad;
    
    qualitymap qm; 
    
    bool showColor;
    
    
};
#endif

