#include "DNA.h"

#include <math.h> 

DNA::DNA(){

}

void DNA::append(string& sDNA){
    base b;
    for(auto it = sDNA.begin(); it != sDNA.end(); ++it) {
        b = char2DNA(*it);
        sequence.push_back(b);
    }
}

void DNA::addQuality(string& qal){
    uint q;
    uint i = 0;
    for(auto it = qal.begin(); it != qal.end(); ++it) {
        q = static_cast<uint>(*it);

        if( i < sequence.size()){
            sequence[i].quality = q;
        }
        i++;
    }
}

base DNA::char2DNA(char& dnaChar){
    base b;
    b.code = dnaChar;

    return b;
}

void DNA::printColoredDNA(WINDOW* win, std::pair<uint, uint> p, bool changeColor){
    base b;
    char a;

    for(auto it = sequence.begin(); it != sequence.end(); ++it) {
        b = (*it);
        a = b.code;

        if(changeColor){
            wattron(win, COLOR_PAIR(b.quality));
            waddch(win, a);
            //wattroff(win, COLOR_PAIR(b.quality));
        }else{
            waddch(win, a);
        }
    }
}



