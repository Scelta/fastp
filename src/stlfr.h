#ifndef STLFR_H
#define STLFR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <unordered_map>

#include "options.h"
#include "read.h"

using namespace std;

typedef unordered_map<string,int> stringmap;
typedef unordered_map<int,int> intmap;

class stlfr{
public:
    stlfr(Options* opt);
    ~stlfr();
    void ReadBarcodeList(string file);
    void findBarcode(Read* r);
    void process(Read* r1, Read* r2 = NULL);
    void addstlfrToName(Read* r, string stlfr);

    static bool test();

private:
    Options* mOptions;
    long barcodeSpace;
};


#endif
