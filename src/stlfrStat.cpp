#include <stdlib.h>
#include "stlfrStat.h"
#include "stats.h"
#include "htmlreporter.h"
#include <memory.h>

#define MAXBCODE 65535

//init for stlfr stats
StlfrStats::StlfrStats(Options* opt){
  mOptions = opt;
  mStlfrValid = 0;
  mStlfrMulti = 0;
  mStlfrMulti10 = 0;
  mStlfrMulti100 = 0;
  mStlfrMulti255 = 0;

  mB1 = 0;
  mB2 = 0;
  mB3 = 0;
  int sbSize = mOptions->stlfr.barcodeSpace; // not include zero
  mStlfrBarcode = new unsigned short int **[sbSize+1];
  memset(mStlfrBarcode, 0, sizeof(unsigned short int)*(sbSize+1));
  for(int i=0; i<=sbSize; i++){
    mStlfrBarcode[i] = new unsigned short int *[sbSize+1];
    memset(mStlfrBarcode[i], 0, sizeof(unsigned short int)*(sbSize+1));
    for(int j=0; j<=sbSize; j++){
      mStlfrBarcode[i][j] = new unsigned short int [sbSize+1];
      memset(mStlfrBarcode[i][j], 0, sizeof(unsigned short int)*(sbSize+1));
    }
  }
}

StlfrStats::~StlfrStats(){
  //delete array
  int sbSize = mOptions->stlfr.barcodeSpace;
  for(int i=0; i<=sbSize; i++){
    for(int j=0; j<=sbSize; i++){
      delete mStlfrBarcode[i][j];
    }
    delete mStlfrBarcode[i];
  }
  delete mStlfrBarcode;
}

void StlfrStats::statStlfr(Read* r){
  int b1 = r->stlfrB1;
  int b2 = r->stlfrB2;
  int b3 = r->stlfrB3;
  if( b1 > mB1 )
    mB1 = b1;
  if( b2 > mB2 )
    mB2 = b2;
  if( b3 > mB3 )
    mB3 = b3;

  if(mStlfrBarcode[b1][b2][b3] < MAXBCODE){
    mStlfrBarcode[b1][b2][b3] ++;
  }
}

//merge stlfr barcode frequency
StlfrStats* StlfrStats::merge(vector<StlfrStats*>& list) {
  //cout << "list.size()： " << list.size() << endl;
  if(list.size() == 0)
      return NULL;

  // barcode stats
  StlfrStats* s = new StlfrStats(list[0]->mOptions);

  for(int t=0; t<list.size(); t++) {
    if(list[t]->mB1 > s->mB1)
        s->mB1 = list[t]->mB1;
    if(list[t]->mB2 > s->mB2)
        s->mB2 = list[t]->mB2;
    if(list[t]->mB3 > s->mB3)
        s->mB3 = list[t]->mB3;

    for(int b1=1; b1<=s->mB1; b1++) {
      for(int b2=1; b2<=s->mB2; b2++) {
        for(int b3=1; b3<=s->mB3; b3++) {

          if(list[t]->mStlfrBarcode[b1][b2][b3] > 0){
            unsigned int sCount = s->mStlfrBarcode[b1][b2][b3];
            unsigned int tCount = list[t]->mStlfrBarcode[b1][b2][b3];
            if(sCount + tCount < MAXBCODE){
              s->mStlfrBarcode[b1][b2][b3] += list[t]->mStlfrBarcode[b1][b2][b3];
            }else {
              s->mStlfrBarcode[b1][b2][b3] = MAXBCODE;
            }
          }
          // stat
          if(t == (list.size() - 1)){
            if(s->mStlfrBarcode[b1][b2][b3] > 0){
              s->mStlfrValid ++;
              if(s->mStlfrBarcode[b1][b2][b3] > 1){
                s->mStlfrMulti ++;
                if(s->mStlfrBarcode[b1][b2][b3] >= 10){
                  s->mStlfrMulti10 ++;
                  if(s->mStlfrBarcode[b1][b2][b3] >= 100){
                    s->mStlfrMulti100 ++;
                    if(s->mStlfrBarcode[b1][b2][b3] >= MAXBCODE){
                      s->mStlfrMulti255 ++;
                    }
                  }
                }
              }
            }
          }
          // end stat
        }
      }
    }
  }
  return s;
}

void StlfrStats::print() {
    cerr << "Valid barcode found: " << mStlfrValid << endl;
    cerr << "Valid & hit >  10^0: " << mStlfrMulti << endl;
    cerr << "Valid & hit >= 10^1: " << mStlfrMulti10 << endl;
    cerr << "Valid & hit >= 10^2: " << mStlfrMulti100 << endl;
    cerr << "Valid & hit >=  MAX: " << mStlfrMulti255 << endl;

}

void StlfrStats::reportJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    // barcode frequency
    ofs << padding << "\t" << "\"stLFR_barcode_frequency\": {" << endl << padding << "\t";
    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    //printf("mB1~mB3: %d,%d,%d\n",mB1,mB2,mB3);
    for(int b1=1; b1<=mB1; b1++) {
      count1 ++;
      bool empty = true;
      for(int b2=1; b2<=mB2; b2++) {
        count2 ++;
        for(int b3=1; b3<=mB3; b3++) {
          count3 ++;
          int count = mStlfrBarcode[b1][b2][b3];
          //printf("count: %d, mStlfrBarcode[b1][b2][b3]: %d\n",count,mStlfrBarcode[b1][b2][b3]);
          if(count > 0){
            empty = false;
            char tag[16];
            sprintf(tag,"\t\"%04d_%04d_%04d\":%d", b1, b2, b3, count);
            ofs << tag;
            if(count1 != mB1 && count2 != mB2 && count3 != mB3)
              ofs << ",";
          }
        }
      }
      if(empty==false)
        ofs << endl << padding << "\t";
    }
    ofs << padding << "\t" << "}," << endl;

    ofs << padding << "}," << endl;
}
