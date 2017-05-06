#include <iostream>
#include <stdio.h>
#include <curses.h>
#include <list>
//#include <string.h>
#include <stdlib.h>
#include <math.h> 
#include <vector>
#include <bitset>
#include <map>
#include <future>

using namespace std;

#include "options.h"
#include "DNA.h"
#include "fastq.h"
#include "fqless.h"

string version = "Version 0.2 Alpha";




void showTheHelp(){
    cout << "fqless is a less like tool for fastq files" << std::endl;
    cout << version << std::endl;
    cout << "" << std::endl;
    cout << "usage: fqless [-b 100] [filename]"  << std::endl;
    cout << "" << std::endl;
    cout << "-b, --buffer   integer, specifying number of screens to keep in" << std::endl;
    cout << "               buffer for scrolling (default: 100, min: 5)" << std::endl;
    cout << "-d, --debug    enable some debugging info" << std::endl;
    cout << "Color support:" << std::endl;
    cout << "The main feature of fqless is to hide the quality line of a fastq file" << std::endl;
    cout << "and instead color code the DNA accordingly." << std::endl;
    cout << "Activate 256-color support for your terminal to use this." << std::endl;
}


void showTheVersion(){
    cout << version << std::endl;
}

int main(int argc, char * argv[]) {

    options * opts    = new options();
    opts->qm          = buildQualityMap();
    opts->offset      = 0;
    opts->tellg       = 0;
    opts->firstInPad  = 0;
    opts->lastInPad   = -1;
    opts->qualitycode = 0;
    opts->buffersize  = 100; // must be at least 3
    opts->showColor   = true;
    opts->debug       = false;


    // lets define options and defaults
    char* input         = NULL;
    int c;

    int showHelp    = 0;
    int showVersion = 0;
    string outputDir;


    const struct option longopts[] = {
        {"version",   no_argument,        0, 'v'},
        {"help",      no_argument,        0, 'h'},
        {"debug",      no_argument,        0, 'd'},
        {"buffer", required_argument,     0, 'b'},
        {0,0,0,0},                          
    };

    int option_index = 0;

    //turn off getopt error message
    opterr=0; 


    while((c = getopt_long (argc, argv, ":hvdb:W",  longopts, &option_index)) != -1){ 
        switch (c) {
            case 'b':
                opts->buffersize = atoi(optarg);
                break;
            case 'h':
                showHelp = 1;
                break;
            case 'v':
                showVersion = 1;
                break;
            case 'd':
                opts->debug = true;
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
    // input has to be the last argument
    if(argc > 1){
        input = argv[argc-1];
    }


    // minimal value for buffersieze is 5
    if(opts->buffersize < 5){
        opts->buffersize = 5;
    }

    if(showHelp){
        showTheHelp();
        exit(0);
    }
    if(showVersion){
        showTheVersion();
        exit(0);
    }
    
    

    if(input == NULL){
        cerr << "No input given. See fqless -h for help." << std::endl;
        exit(0);
    }else{
        std::ifstream infile(input);
        if(infile.good() != true){
            cerr << "Can not find or open file: " << input << std::endl;
            exit(0);
        }
    }
    opts->input = input;

    // init a new fql less class
    fqless fql = fqless(opts);

    // exit
    exit(0);


}

