#include "DNA.h"

/*
dnamap DNA::getDNAMap(){	
	map <char, Sdna> dnaMap;	
	// A is bitflip T and G can be flipped to C
	dnaMap['A'] = bitset<3>(001);
	dnaMap['T'] = bitset<3>(110);
	dnaMap['G'] = bitset<3>(011);
	dnaMap['C'] = bitset<3>(100);
	
//	dnaMap['N'] = bitset<3>(111);
	dnaMap['N'] = bitset<3>(000);
	return dnaMap; 
}
*/

/*
bitmap DNA::getBitMap(){
	bitmap revMap;
	cout << "DNA Bitmap";
	//for(map<char, Sdna>::iterator i = dnaMap.begin(); i != dnaMap.end(); i++){
	//	revMap[*(i)->second] = *(i)->first;
	//	cout << i->first;
	//}
	//for (auto i = dnaMap.begin(); i != dnaMap.end(); ++i)
     //   revMap[i->second] = i->first;

	return revMap;
}
*/
/*
dnamap DNA::getDNAMap(){
	dnamap dnaMap;
	baseT baseA = 'A';
	baseT baseC = 'C';
	baseT baseG = 'G';
	baseT baseT = 'T';
	
	
	dnaMap.insert(Sdna(001), baseA);
	dnaMap.insert(Sdna(011), baseC);
	dnaMap.insert(Sdna(101), baseG);
	dnaMap.insert(Sdna(111), baseT);
//	m.insert(dnamap::value_type("A", bitset<3>(001)));
	//bob.insert("A",1);
	//return bob;
	return dnaMap;
}*/
/*
DNA::DNA(vector <char> dnaVec){
    int len = dnaVec.size();
//    if(len == 0){
//        cerr << "DNA of length 0"; 
        //std::exit(0);
        
  //  }
    vector < char > test;
    //return  test;
}*/

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


string DNA::getSequence(){
	string DNAstring;
	uint i = 0;
	
	for(auto it = sequence.begin(); it != sequence.end(); it++){
		DNAstring[i] = DNA2char(*it);
		i++;
	}	

	return DNAstring;
}
/*
Tdna DNA::getSequence(){
	return sequence;
}*/




base DNA::char2DNA(char& dnaChar){
    base b;
    switch(dnaChar) {
        case 'A' : 
                 b.code = Sdna(1);
                 break;
        case 'C' :
                 b.code = Sdna(10);
                 break;
        case 'G' :
                 b.code = Sdna(11);
                 break;
        case 'T' : 
                 b.code = Sdna(110);
                 break;
        case 'N' :
                 b.code = Sdna(0);
                 break;
        default :
                 cout << "whats that for a letter?";
                 break;
    }
    
    
	return b;
}
char DNA::DNA2char(base& base){
	char a;
	Sdna b;
	b = base.code;
	
	Sdna A = Sdna(1);
	Sdna C = Sdna(10);
	Sdna G = Sdna(11);
	Sdna T = Sdna(110);
	Sdna N = Sdna(0);
	
	if( b == A ){
	    a = 'A';
	}else if( b == G ){
	    a = 'G';
	}else if( b == C ){
	    a = 'C';
	}else if( b == T ){
	    a = 'T';
	}else if( b == N ) {
	    a = 'N';
	}

	return a;	
}


void DNA::printColoredDNA(WINDOW* win){
    base b;
    char a;
    
    bool changeColor;
    changeColor = can_change_color();

	for(auto it = sequence.begin(); it != sequence.end(); ++it) {
		b = (*it);
		a = DNA2char(b);
		
		if(changeColor){
		    //printw("yes");
		    color c;
		    c = IntToColor(b.quality);
		    //stringstream convert;
		    //convert << "qul:" << b.quality << " " << c.R << " " << c.G << " " <<  c.B << "\n" ;


		    init_color(b.quality, c.R,c.G,c.B);
		    init_pair(b.quality, b.quality, -1);
		    wattron(win, COLOR_PAIR(b.quality));
		    waddch(win, a);
            wattroff(win, COLOR_PAIR(b.quality));
		}else{
		    printw("no");
		}
	}
}



color DNA::IntToColor(int i){
        // these are the maximal ASCII values used for quality scores
        // we can map everything into this range
        uint max;
        uint min;
        float range;
        max     = 126;
        min     = 33;
        range   = 41;
        
        
        // floor the number, so 0 is 0.
        float n;
        n       = i - min;
        n       /= range;
        max     = max - min;
        
       
        color result;
        float m  = 1;
        uint scale = 1000;
        result.R = floor(( (1 - n) * m ) * scale);
        result.G = floor(( n * m )* scale);
        result.B = 0;
        
        if(result.R > scale){ result.R = scale;}
        if(result.G > scale){ result.G = scale;}
        
        //cout << n << " " << result.R << " " << result.G << "- " ;
     
        return result;
    }
