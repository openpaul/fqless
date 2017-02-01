#include <stdio.h>
#include "options.h"




//#include "DNA.h"
#include "fastq.h"

string version = "Version 0.1 Alpha";


WINDOW * Wtext;
WINDOW * Wcmd;
WINDOW * lNumb;






int buffersize = 2;

void quit(){
    delwin(Wtext);
    delwin(Wcmd);
    endwin();
}

void winInit(options * opts){

    opts->textrows = LINES-5;
    if(opts->linenumbers){     
        opts->textcols = COLS-1;
    }else{
        opts->textcols = COLS;
    }
    Wcmd            = newwin(1, COLS, opts->textrows+4, 0);
    Wtext           = newpad(buffersize*opts->textrows, opts->textcols);
    
    // init line number pad
    if(opts->linenumbers){
        lNumb       = newpad(buffersize*opts->textrows, opts->linenumberspace);
    }
    
    keypad(Wcmd, TRUE);
}



void fillPad(options* opts, int offset, fastq* FQ, int dir=1){
    int b;
    int bmore;
    opts->avaiLines   = buffersize*opts->textrows;
    bmore               = FQ->readmore(offset, opts->avaiLines, opts->textcols, dir, Wtext);
    
    // we may need to make space for a sequence that a little bit larger then 
    // the buffer and thus then we expected
    wprintw(Wtext, "Lines %i ", opts->avaiLines);
    //wprintw(Wtext, "Lines %i and we req %i, so resize\n", bmore, opts->avaiLines);
    if(bmore != opts->avaiLines and bmore > 0){
      // wprintw(Wtext, "We should ");
        opts->avaiLines = bmore;
        wresize(Wtext, opts->avaiLines, opts->textcols);
        wresize(lNumb, opts->avaiLines, opts->linenumberspace);
       // wprintw(Wtext, "COLS %i ROWS %i\n", opts->linenumberspace, opts->avaiLines);

        //opts->textrows = LINES-5;
        wprintw(Wtext, "Lines %i ", opts->avaiLines);
    }


    int i = 0; // keep track of lines
    for(auto it = FQ->content.begin(); it != FQ->content.end(); ++it) {
        fastqSeq fq = *it;
        
        if(i > 0) wprintw(Wtext, "\n"); // spacer
        // paint the name
        wattron(Wtext, COLOR_PAIR(1));
        wprintw(Wtext, fq.name.c_str());
	    wattroff(Wtext, COLOR_PAIR(1));
        wprintw(Wtext, "\n");

        // print sequence
        fq.dna.printColoredDNA(Wtext);
        wprintw(Wtext, "\n");
        
        
        i = i + 2 + ceil(fq.dna.sequence.size()/COLS) ;
    }

     
      //wprintw(Wtext, "AHHHH %i ", opts->avaiLines);
     // set line numbers
     for(uint i = 1; i <= opts->avaiLines; i++){
        wattron(lNumb, COLOR_PAIR(1));
        wprintw(lNumb, "%i\n", i);
     }
     wattroff(lNumb, COLOR_PAIR(1));
     
}
 



int main(int argc, char * argv[]) {
    
    options* opts;
    opts->linenumbers = true;
    if(opts->linenumbers){
        opts->linenumberspace = 5;
    }else{
        opts->linenumberspace = 0;
    }
    
    opts->offset = 0;
    
    

    initscr();                      // start ncurses
    curs_set(0);                    // hide cursor
    start_color();                  // start color mode
    assume_default_colors(-1,-1);   // make transparent mode possible if supported
    atexit(quit);                   // what to do at exit
    init_pair(1, -1, -1);           // default colors, white on standart background

    winInit(opts);
    
    // keep track of offset of file
    int offset;
    offset = 0;

    
    // lets define options and defaults
    char* input;
    int do_verbose, do_help, c;
    string outputDir;

    
    
    
    const struct option longopts[] =
    {
        {"version",   no_argument,        0, 'v'},
        {"help",      no_argument,        0, 'h'},
        {"input",     required_argument,  NULL, 'i'},
        {0,0,0,0},
    };

    int index;
    int iarg=0;

    //turn off getopt error message
    opterr=1; 

    while((c = getopt_long(argc, argv, ":hvi::W;", longopts, NULL)) != -1){ 
        switch (c) {
        case 'i':
            input = optarg;
            break;
        case 'h':
            do_help = 1;
        case 'v':
            do_verbose = 1;
            break;
        case 0:     
            break;

        case 1:
        break;
        case ':': 
            fprintf(stderr, "%s: option `-%c' requires an argument\n",
                    argv[0], optopt);
            break;
        case '?':
        default:  
            fprintf(stderr, "%s: option `-%c' is invalid: ignored\n",
                    argv[0], optopt);
            break;
        }
    }


    
    if(input != NULL){
        fastq* FQ = new fastq(input);
        
        // update status
        mvwprintw(Wcmd, 0,0, "File contains %d Sequences ", FQ->nOfSequences);

        
        // update pad
        fillPad(opts, 0, FQ);
        
        refresh();
        prefresh(Wtext, offset,0,0, opts->linenumberspace, opts->textrows, opts->textcols);
        prefresh(lNumb, offset,0,0, 0, opts->textrows, opts->linenumberspace);
        wrefresh(Wcmd);

        int ch;

        noecho();
        while(1) {
            ch = wgetch(Wcmd);
            

            //wechochar(Wtext, ch);
            //wrefresh(Wtext);
            switch(ch){
            case KEY_UP:
	            offset--;
	            if(offset < 0) offset = 0;
	            //refresh();
	            mvwprintw(Wcmd, 0,0, "offset %i lines %i ", offset, buffersize*(LINES-1));
	            pnoutrefresh(Wtext, offset,0,0, opts->linenumberspace, opts->textrows, opts->textcols);
                prefresh(lNumb, offset,0,0, 0, opts->textrows, opts->linenumberspace);
	            break;

            case KEY_DOWN:
	            
	            //refresh();
	            if( offset < (opts->avaiLines - opts->textrows) -1){ // TODO why the 1? lets find out later
    	            offset++;
	                mvwprintw(Wcmd, 0,0, "offset %i lines %i ", offset, buffersize*(LINES-1));
                    pnoutrefresh(Wtext, offset,0,0, opts->linenumberspace, opts->textrows, opts->textcols);  
                    prefresh(lNumb, offset,0,0, 0, opts->textrows, opts->linenumberspace);
                    
                }    
                 
	            break;

            case KEY_NPAGE:

	            break;

            case KEY_PPAGE:
	            //form_driver(form, REQ_PREV_PAGE);
	            //set_form_page(form, --cur_page);
	            break;
	        case KEY_RESIZE:
	            
                endwin();
                winInit(opts);
                fillPad(opts, offset, FQ);
                refresh();
                prefresh(Wtext, offset,0,0, 0, LINES-1, COLS);  
	            break;
            }
        }  
          
        getch();
        //cout << "COLORS" << COLORS;
    } 
}

