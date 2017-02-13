

#include "fastq.h"
#include "fqless.h"
#include "options.h"

string version = "Version 0.1 Alpha";



void fqless::quit(){
    endwin();
}

void showTheHelp(){
    cout << "This would be the help" << std::endl;
}


void showTheVersion(){
    cout << version << std::endl;
}

int main(int argc, char * argv[]) {

    options * opts  =  new options();

    opts->qm   = buildQualityMap();


    opts->offset        = 0;
    opts->tellg         = 0;
    opts->firstInPad    = 0;
    opts->lastInPad     = -1;
    opts->qualitycode   = 0;
    opts->buffersize    = 100;
    opts->showColor     = true;





    // lets define options and defaults
    char* input         = NULL;
    int c, ch;

    int showHelp       = 0;
    int showVersion    = 0;
    string outputDir;


    const struct option longopts[] = {
        {"version",   no_argument,        0, 'v'},
        {"help",      no_argument,        0, 'h'},
        {"input",     required_argument,  NULL, 'i'},
        {0,0,0,0},
    };

    int option_index = 0;

    //turn off getopt error message
    opterr=0; 

    while((c = getopt_long (argc, argv, ":hvi:W",  longopts, &option_index)) != -1){ 
        switch (c) {
            case 'i':
                input = optarg;
                break;
            case 'h':
                showHelp = 1;
                break;
            case 'v':
                showVersion = 1;
                break;
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


    if(showHelp){
        showTheHelp();
        exit(0);
    }
    if(showVersion){
        showTheVersion();
        exit(0);
    }

    opts->input = input;

    // init a new fql less class
    fqless fql = fqless(opts);


    //getch();


}

