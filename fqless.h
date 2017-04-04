
#ifndef FQLESS_H
#define FQLESS_H

class fqless{
public:
    int ch;
    string colorMessage = "";
    WINDOW * Wtext;
    WINDOW * Wcmd;
    // init the class
    fqless(options*); 

    static void quit();
    void initTheColors(std::pair<uint, uint> p);
    void winInit(options * opts);
    color IntToColor(int i, std::pair<uint, uint> p);



private:
    void fillPad(options* opts, fastq* FQ, int dir);
    void statusline(options* , fastq* ,WINDOW* );
};

#endif
