
#ifndef OPTIONS_H
#define OPTIONS_H


#include <string>
#include <map>
#include <vector>

using namespace std;
typedef unsigned int  uint;
typedef std::map<std::string, std::pair<int, int> > qualitymap;





struct options{
    char* input;
    
    bool linenumbers;
    int buffersize; 
   
    
    bool debug;
    
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
    bool basicColor; // use default terminal colors
    
};
#endif

