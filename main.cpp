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

string version = "Version 0.1 Alpha";




void showTheHelp(){
    cout << "fqless is a less like tool for fastq files" << std::endl;
    cout << version << std::endl;
    cout << "" << std::endl;
    cout << "usage: fqless -i fastq.file" << std::endl;
    cout << "" << std::endl;
    cout << "Color support:" << std::endl;
    cout << "A main feature of fqless is to hide the quality line of a fastq file" << std::endl;
    cout << "and instead color code the DNA accordingly." << std::endl;
    cout << "Activate color support for your terminal to use this." << std::endl;
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


    // lets define options and defaults
    char* input         = NULL;
    int c, ch;

    int showHelp    = 0;
    int showVersion = 0;
    int  buffer  = 0;
    string outputDir;


    const struct option longopts[] = {
        {"version",   no_argument,        0, 'v'},
        {"help",      no_argument,        0, 'h'},
        {"buffer", required_argument,     0, 'b'},
        {"input",     required_argument,  NULL, 'i'},
        {0,0,0,0},                          
    };

    int option_index = 0;

    //turn off getopt error message
    opterr=0; 

    while((c = getopt_long (argc, argv, ":hvb:i:W",  longopts, &option_index)) != -1){ 
        switch (c) {
            case 'i':
                input = optarg;
                break;
            case 'b':
                opts->buffersize = atoi(optarg);
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
        cout << "No input given. See fqless -h for help." << std::endl;
        exit(0);
    }

    opts->input = input;

    // init a new fql less class
    fqless fql = fqless(opts);




}

