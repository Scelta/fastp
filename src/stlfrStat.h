#ifndef STLFR_STAT_H
#define STLFR_STAT_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <math.h>
#include <unordered_map>

#include "options.h"
#include "read.h"

using namespace std;

class StlfrStats{
public:
    StlfrStats(Options* opt);
    ~StlfrStats();
    void ReadBarcodeList(string file);
    void findBarcode(Read* r);
    void process(Read* r1, Read* r2 = NULL);
    void addstlfrToName(Read* r, string stlfr);
    void reportJson(ofstream& ofs, string padding);
    void statStlfr(Read* r);
    static StlfrStats* merge(vector<StlfrStats*>& list);
    void print();
    static bool test();

private:
    Options* mOptions;
    unsigned short int ***mStlfrBarcode;
    int mStlfrValid;
    int mStlfrMulti;
    int mStlfrMulti10;
    int mStlfrMulti100;
    int mStlfrMulti255;

    int barcodeSpace;
    int mB1;
    int mB2;
    int mB3;
};


#endif
