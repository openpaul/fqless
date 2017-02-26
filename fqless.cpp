#include "DNA.h"
#include "fastq.h"
#include "fqless.h"
#include <math.h>
#include <future>

color fqless::IntToColor(int i, std::pair<uint, uint> p){
    // these are the maximal ASCII values used for quality scores
    // we can map everything into this range
    uint max;
    uint min;
    float range;

    min   = p.first;
    max   = p.second;
    range = max - min;

    float n;
    n       = i - min;
    n      /= range;


    color result;
    float m     = 1;
    uint scale  = 1000;
    result.R    = floor(( (1 - n) * m ) * scale);
    result.G    = floor(( n * m )* scale);
    result.B    = 0;

    if(result.R > scale){ result.R = scale;}
    if(result.G > scale){ result.G = scale;}

    return result;
}

void fqless::initTheColors(std::pair<uint, uint> p){
    int fqmin, fqmax;
    fqmin = p.first;
    fqmax = p.second;
    color c;
    while(fqmin < fqmax){
        c = IntToColor(fqmin, p);
        init_color(fqmin, c.R,c.G,c.B);
        init_pair(fqmin, fqmin, -1);
        
        fqmin++;
    }
}

void fqless::winInit(options * opts){

    opts->textrows = LINES-2;
    opts->textcols = COLS;

    Wcmd            = newwin(1, COLS, opts->textrows+1, 0);
    Wtext           = newpad(opts->buffersize*opts->textrows, opts->textcols);



    keypad(Wcmd, TRUE);
}


void fqless::quit(){
    endwin();
}

void fqless::fillPad(options* opts, fastq* FQ, int dir=1){

    // clean the pad now
    wclear(Wtext); 

    std::pair<uint, uint> qalpair;
    if(opts->showColor){
        qalpair = opts->qm.at(FQ->possibleQual[opts->qualitycode]);
    }

    // we may need to make space for a sequence that a little bit larger then 
    // the buffer and thus then we expected

    //wprintw(Wtext, "Index size %i \n",  FQ->index.size());
    if(opts->linesTohave != opts->avaiLines){
        opts->avaiLines = opts->linesTohave;
        wresize(Wtext, opts->avaiLines+1, opts->textcols);
    }


    int i = 0; // keep track of lines
    for(auto& it: FQ->content) {
        fastqSeq& fq = it;


        if(i > 0) wprintw(Wtext, "\n"); // spacer
        // paint the name
        wattron(Wtext, COLOR_PAIR(1));
        wprintw(Wtext, fq.name.c_str());
        wattroff(Wtext, COLOR_PAIR(1));
        wprintw(Wtext, "\n");

        // print sequence

        fq.dna.printColoredDNA(Wtext, qalpair, opts->showColor);
        wprintw(Wtext, "\n");
        //fq.inpad    = true;

        i = i + 2 + ceil(fq.dna.sequence.size()/COLS) ;
    }


}

fqless::fqless(options* opts){

    if(opts->input != NULL){

        initscr();                      // start ncurses
        curs_set(0);                    // hide cursor
        start_color();                  // start color mode
        assume_default_colors(-1,-1);   // make transparent mode possible if supported
        atexit(quit);                   // what to do at exit
        init_pair(1, -1, -1);           // default colors, white on standart background

        winInit(opts);

        noecho();             
        opts->avaiLines   = opts->buffersize*opts->textrows;

        // switch to black and white if we can not show color
        if(can_change_color() == false){
            opts->showColor = false;
        }

        fastq* FQ = new fastq(opts->input);

        // here we can build the first index
        // it should index as much as the first buffer only
        FQ->buildIndex(opts);
        FQ->showthese(opts, 1, Wtext);

        opts->offset        = 0;

        if(opts->showColor == true){
            initTheColors(opts->qm.at(FQ->possibleQual[opts->qualitycode]));  
        }else{
            if(can_change_color() == false){
                colorMessage = " | no color support by the terminal";
            }else{
                colorMessage = " | no valid quality range found (" + to_string(FQ->minQal) + ","  + to_string(FQ->maxQal) + ")";
            }
        }




        // update pad
        fillPad(opts, FQ);

        refresh();
        pnoutrefresh(Wtext, opts->offset,0,0, 0, opts->textrows, opts->textcols);
        wrefresh(Wcmd);


        // launcyh an async thread to build an index
        // async because, we do not need this to continue 
        // but it is needed later to quickly load more data
        auto indexThread = std::async(std::launch::async, &fastq::buildIndex, FQ, opts);


        // update status
        mvwprintw(Wcmd, 0,0, "File: %s%s", FQ->file.c_str(),colorMessage.c_str());

        while(1) {
            ch = wgetch(Wcmd);

            switch(ch){
                case KEY_UP:
                    opts->offset--;
                    if(opts->offset < 0 and opts->firstInPad > 0) {
                        FQ->showthese(opts, 0, Wtext);
                        opts->offset--;
                        fillPad(opts, FQ);
                    }
                    if(opts->offset < 0){opts->offset = 0;} 


                    //mvwprintw(Wcmd, 0,0, "offset %i lines %i possible offs %i ", opts->offset, opts->buffersize*opts->textrows, opts->avaiLines - opts->textrows -1 );
                    pnoutrefresh(Wtext, opts->offset,0,0, 0, opts->textrows, opts->textcols);
                    doupdate();
                    break;

                case KEY_DOWN:
                    opts->offset++;

                    if( (int)opts->offset > (int) (opts->avaiLines - opts->textrows -2) && ((int)  opts->lastInPad < (int) FQ->index.size())){
                        FQ->showthese(opts, 1, Wtext);
                        fillPad(opts, FQ);
                    }
                    if((int)opts->offset > (int)(opts->avaiLines - opts->textrows - 1)){
                        opts->offset = opts->avaiLines - opts->textrows -1;
                    }
                    //mvwprintw(Wcmd, 0,0, "offset %i lines %i possible offs %i ", opts->offset, opts->buffersize*(LINES-1), opts->avaiLines - opts->textrows -1 );
                    pnoutrefresh(Wtext, opts->offset,0,0, 0, opts->textrows, opts->textcols);                          
                    doupdate();
                    break;

                case KEY_RIGHT:
                    if(opts->showColor == false){
                        break;
                    }
                    opts->qualitycode++;
                    if(opts->qualitycode > (int)FQ->possibleQual.size() -1 ){
                        opts->qualitycode = 0;
                    }
                    // colors changed
                    initTheColors(opts->qm.at(FQ->possibleQual[opts->qualitycode]));
                    fillPad(opts, FQ);
                    pnoutrefresh(Wtext, opts->offset,0,0, 0, opts->textrows, opts->textcols);  
                    doupdate();
                    break;

                case KEY_LEFT:
                    if(opts->showColor == false){
                        break;
                    }
                    opts->qualitycode--;
                    if(opts->qualitycode < 0){
                        opts->qualitycode = FQ->possibleQual.size() - 1;
                    }
                    // colors changed
                    initTheColors(opts->qm.at(FQ->possibleQual[opts->qualitycode]));  
                    fillPad(opts, FQ);
                    pnoutrefresh(Wtext, opts->offset,0,0, 0, opts->textrows, opts->textcols);  
                    doupdate();
                    break;

                case KEY_NPAGE:
                    opts->offset += opts->textrows;

                    if( (int)opts->offset > (int) (opts->avaiLines - opts->textrows -2) && ((int)  opts->lastInPad < (int) FQ->index.size())){
                        FQ->showthese(opts, 1, Wtext);
                        fillPad(opts, FQ);
                    }
                    if((int)opts->offset > (int)(opts->avaiLines - opts->textrows - 1)){
                        opts->offset = opts->avaiLines - opts->textrows -1;
                    }
                    pnoutrefresh(Wtext, opts->offset,0,0, 0, opts->textrows, opts->textcols);  
                    doupdate(); 
                    break;

                case KEY_PPAGE:
                    opts->offset -= opts->textrows;
                    if(opts->offset < 0 and opts->firstInPad > 0) {
                        FQ->showthese(opts, 0, Wtext);
                        opts->offset--;
                        fillPad(opts, FQ);
                    }
                    if(opts->offset < 0){opts->offset = 0;} 

                    pnoutrefresh(Wtext, opts->offset,0,0, 0, opts->textrows, opts->textcols);
                    doupdate();
                    break;

                case KEY_RESIZE:

                    endwin();
                    winInit(opts);
                    fillPad(opts, FQ);
                    refresh();
                    prefresh(Wtext, opts->offset,0,0, 0, LINES-1, COLS);  
                    break;

                case 'q':
                    quit();
                    exit(0);
                    break;
            }

        }  

    }
    return;
}

