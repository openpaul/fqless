#include <stdio.h>





//#include "DNA.h"
#include "fastq.h"

string version = "Version 0.1 Alpha";


WINDOW * Wtext;
WINDOW * Wcmd;
WINDOW * lNumb;






int buffersize = 50;

void quit(){
    //delwin(Wtext);
    //delwin(Wcmd);
    endwin();
}

void showTheHelp(){
    cout << "This would be the help" << std::endl;
}


void showTheVersion(){
    cout << version << std::endl;
}

void winInit(options * opts){

    opts->textrows = LINES-2;
    if(opts->linenumbers){     
        opts->textcols = COLS-5;
    }else{
        opts->textcols = COLS;
    }
    Wcmd            = newwin(1, COLS, opts->textrows+1, 0);
    Wtext           = newpad(buffersize*opts->textrows, opts->textcols);
    
    // init line number pad
    if(opts->linenumbers){
        lNumb       = newpad(buffersize*opts->textrows, opts->linenumberspace);
    }
    
    keypad(Wcmd, TRUE);
}



void fillPad(options* opts, fastq* FQ, int dir=1){
    
    // clean the pad now
    wclear(Wtext); 
    wclear(lNumb);
    
    // we may need to make space for a sequence that a little bit larger then 
    // the buffer and thus then we expected

    //wprintw(Wtext, "Index size %i \n",  FQ->index.size());
    if(opts->linesTohave != opts->avaiLines){
        opts->avaiLines = opts->linesTohave;
        wresize(Wtext, opts->avaiLines+1, opts->textcols);
        wresize(lNumb, opts->avaiLines+1, opts->linenumberspace);
    }


    int i = 0; // keep track of lines
    int j = 1; // offset in line numbers
    for(auto& it: FQ->content) {
        fastqSeq& fq = it;
        
        //if(fq.inpad == false){
        //    //fq.inpad = false;            
        //    j = j + 2 + ceil(fq.dna.sequence.size()/COLS) ;
        //    continue;
        //}
                
        if(i > 0) wprintw(Wtext, "\n"); // spacer
        // paint the name
        wattron(Wtext, COLOR_PAIR(1));
        wprintw(Wtext, fq.name.c_str());
	    wattroff(Wtext, COLOR_PAIR(1));
        wprintw(Wtext, "\n");

        // print sequence
        fq.dna.printColoredDNA(Wtext);
        wprintw(Wtext, "\n");
        //fq.inpad    = true;
        
        i = i + 2 + ceil(fq.dna.sequence.size()/COLS) ;
    }

     

     // set line numbers
     for(int i = j; i <= opts->avaiLines; i++){
        wattron(lNumb, COLOR_PAIR(1));
        wprintw(lNumb, "%i\n", i);
     }
     wattroff(lNumb, COLOR_PAIR(1));
       
}
 



int main(int argc, char * argv[]) {
   
    options * opts  =  new options();
    
    opts->linenumbers = false;
    
    if(opts->linenumbers){
        opts->linenumberspace = 5;
    }else{
        opts->linenumberspace = 0;
    }
 
    opts->offset        = 0;
    opts->tellg         = 0;
    opts->firstInPad    = 0;
    opts->lastInPad     = -1;

    
    
     // lets define options and defaults
    char* input;
    int showVersion,showHelp, c, ch;
    string outputDir;

    
    const struct option longopts[] =
    {
        {"version",   no_argument,        0, 'v'},
        {"help",      no_argument,        0, 'h'},
        {"input",     required_argument,  NULL, 'i'},
        {0,0,0,0},
    };

    //int index;
    //int iarg=0;

    //turn off getopt error message
    opterr=1; 

    while((c = getopt_long(argc, argv, ":hvi:W;", longopts, NULL)) != -1){ 
        switch (c) {
        case 'i':
            input = optarg;
            break;
        case 'h':
            showHelp = 1;
        case 'v':
            showVersion = 1;
        case 0:     
            break;
        case 1:
            break;
        case ':': 
            fprintf(stderr, "%s: option `-%c' requires an argument\n",
                    argv[0], optopt);
            break;
        default:  
            fprintf(stderr, "%s: option `-%c' is invalid: ignored\n",
                    argv[0], optopt);
            break;
        }
    }

    
    if(showHelp == 1){
       // quit();
        showTheHelp();
        exit(0);
    }
    if(showVersion == 1){
       // quit();
        showTheVersion();
        exit(0);
    }

   
    


    
   
    
    if(input != NULL){
    
        initscr();                      // start ncurses
        curs_set(0);                    // hide cursor
        start_color();                  // start color mode
        assume_default_colors(-1,-1);   // make transparent mode possible if supported
        atexit(quit);                   // what to do at exit
        init_pair(1, -1, -1);           // default colors, white on standart background
        
        winInit(opts);  
        noecho();             
        opts->avaiLines   = buffersize*opts->textrows;
    
    
        fastq* FQ = new fastq(input);
        // here we can build the first index
        // it should index as much as the first buffer only
        FQ->buildIndex(opts);
        FQ->showthese(opts, 1, Wtext);
        opts->offset        = 0;
        // then we read the first buffer and show it
        
        
        
        // update status
        mvwprintw(Wcmd, 0,0, "File contains %d Sequences ", FQ->nOfSequences);

        
        // update pad
        fillPad(opts, FQ);
        
        
        refresh();
        prefresh(Wtext, opts->offset,0,0, opts->linenumberspace, opts->textrows, opts->textcols);
        prefresh(lNumb, opts->offset,0,0, 0, opts->textrows, opts->linenumberspace);
        wrefresh(Wcmd);
        
        
        // launcyh an async thread to build an index
        // async because, we do not need this to continue 
        // but it is needed later to quickly load more data
        auto indexThread = std::async(std::launch::async, &fastq::buildIndex, FQ, opts);


        

        
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
	                
	                mvwprintw(Wcmd, 0,0, "offset %i lines %i ", opts->offset, buffersize*(LINES-1));
	                pnoutrefresh(Wtext, opts->offset,0,0, opts->linenumberspace, opts->textrows, opts->textcols);
                    prefresh(lNumb, opts->offset,0,0, 0, opts->textrows, opts->linenumberspace);
	                break;

                case KEY_DOWN:
	                opts->offset++;

	                if( opts->offset > (opts->avaiLines - opts->textrows -2) && (opts->lastInPad < FQ->index.size())){
	                    FQ->showthese(opts, 1, Wtext);
	                    fillPad(opts, FQ);
	                }
	                if(opts->offset > (opts->avaiLines - opts->textrows - 1)){
	                    opts->offset = opts->avaiLines - opts->textrows -1;
	                }
	                    mvwprintw(Wcmd, 0,0, "offset %i lines %i possible offs %i ", opts->offset, buffersize*(LINES-1), opts->avaiLines - opts->textrows -1 );
                        pnoutrefresh(Wtext, opts->offset,0,0, opts->linenumberspace, opts->textrows, opts->textcols);  
                        prefresh(lNumb, opts->offset,0,0, 0, opts->textrows, opts->linenumberspace);
                        

                     
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
                    fillPad(opts, FQ);
                    refresh();
                    prefresh(Wtext, opts->offset,0,0, 0, LINES-1, COLS);  
	                break;
            }
        }  
          
        getch();

    } 
}

