#ifndef stlfr_PROCESSOR_H
#define stlfr_PROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include "read.h"

using namespace std;

class stLFRProcessor{
public:
    stLFRProcessor(Options* opt);
    ~stLFRProcessor();
    void process(Read* r1, Read* r2 = NULL);
    void addstLFRToName(Read* r, string stLFR);
    static bool test();

private:
    Options* mOptions;
};


#endif
