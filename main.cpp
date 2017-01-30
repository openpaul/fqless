#include <stdio.h>





//#include "DNA.h"
#include "fastq.h"

string version = "Version 0.1 Alpha";


WINDOW *textfield; // main window
WINDOW *commandline; // make windows possible


void quit(){
    delwin(textfield);
    delwin(commandline);
    endwin();
}

void nextLine(){
    //int x, y;
    // move curser to next line
    //getyx(stdscr, y, x);
    //printw("Koordinatenursprung:       [%d, %d]", y, x);
    //move(y+1, 0);
    

}
 



int main(int argc, char * argv[]) {
    
   // int x, y;
    

    initscr();
    start_color();
//    use_default_colors();
    assume_default_colors(-1,-1);
      atexit(quit);
    //init_color(COLOR_WHITE, 1000,1000,1000);
   //init_pair(1, COLOR_WHITE, COLOR_BLACK);
    init_pair(1, COLOR_WHITE, -1);
  //  keypad(stdscr, TRUE);
    // draw welcome display
   // curs_set(0);
   
   // set sscrolling region

   
 
    
    
    commandline = newwin(1, COLS, LINES-1, 0);
    textfield   = newwin(LINES-1, COLS, 0, 0);
    wsetscrreg(textfield, 0, LINES-1);
    scrollok(textfield, TRUE);
//    
   // mvwprintw(textfield,3, 5, "LINES: %d", LINES);
   // mvprintw(4, 5, "COLS:  %d", COLS);

    //bkgd(COLOR_PAIR(1));
    //wbkgd(textfield, COLOR_PAIR(2));
    //refresh();
     //  wrefresh(commandline);
      //  wrefresh(textfield);;
/* 
initscr();			
 
	if(has_colors() == FALSE)
	{	endwin();
		printf("Your terminal does not support color\n");
		exit(1);
	}
	start_color();			


    init_color(COLOR_RED, 100, 100, 0);
    
	init_pair(1, COLOR_RED, COLOR_BLACK);

	attron(COLOR_PAIR(1));
	print_in_middle(stdscr, LINES / 2, 0, 0, "Viola !!! In color ...");
	attroff(COLOR_PAIR(1));
    	getch();
	endwin();




    getyx(stdscr, y, x);
    mvprintw(5, 5, "Momentane Cursorposition:  [%d, %d]", y, x);
    printw("A");
    printw("A");
    printw("A");
    printw("A");
    printw("A");
    printw("A");
    printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");printw("A");

    getbegyx(stdscr, y, x);
    mvprintw(6, 5, "Koordinatenursprung:       [%d, %d]", y, x);

    getmaxyx(stdscr, y, x);
    mvprintw(7, 5, "Fenstergrße:              [%d, %d]", y, x);

    mvaddstr(11, 2, "Taste drücken -> Ende");

    



    

    */
    
   // cout << "fastq reader " << version << std::endl;
   // cout << std::endl;
    
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
        mvwprintw(commandline, 0,0, "File contains %d Sequences ", FQ->nOfSequences);
/*
        int ch;
        while ((ch = getch()) != KEY_F(1)) {
            switch (ch) {
            case KEY_DOWN:
	            //form_driver(form, REQ_NEXT_FIELD);
	            //form_driver(form, REQ_END_LINE);
	            break;

            case KEY_UP:
	            //form_driver(form, REQ_PREV_FIELD);
	            //form_driver(form, REQ_END_LINE);
	            break;

 
            case KEY_NPAGE:
	            //form_driver(form, REQ_NEXT_PAGE);
	            //set_form_page(form, ++cur_page);
	            break;

            case KEY_PPAGE:
	            //form_driver(form, REQ_PREV_PAGE);
	            //set_form_page(form, --cur_page);
	            break;
            }
        }*/
        
        
        for(uint i = 0; i < FQ->content.size(); i++){
            // print name

		    wattron(textfield, COLOR_PAIR(1));
            const string name = FQ->content[i].name;
		    wattron(textfield, COLOR_PAIR(1));
            wprintw(textfield, "\n");
            wprintw(textfield, name.c_str());
            wprintw(textfield, "\n");

            // print sequence
            FQ->content[i].dna.printColoredDNA(textfield);
            wprintw(textfield, "\n");
        }
        refresh();    
        wrefresh(commandline);
        wrefresh(textfield);
          
        getch();
        //cout << "COLORS" << COLORS;
    } 
}

