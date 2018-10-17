#include "stlfr.h"

stringmap sfLFRbarcodeMap, sfLFRbarcodeSnpMap;
intmap missScore, stlfrLog;
char stlfrBarcodeTag[16];

stlfr::stlfr(Options* opt){
    mOptions = opt;
}

stlfr::~stlfr(){
}

void stlfr::ReadBarcodeList(string file)
{
  //stringmap sfLFRbarcodeMap, sfLFRbarcodeSnpMap;
  //cout << "Run into stlfr::ReadBarcodeList. file: "<< file << endl;
  ifstream fin(file);
  string code;
  int id;
  string base4[4] = {"A","T","C","G"};

  while( fin >> code )
  {
    fin >> id;
    sfLFRbarcodeMap[code] = id;
    //create 1 snp sequence map
    for(int i=0;i<code.length();i++){
      for(int b=0;b<4;b++){
        string snp = code;
        snp.replace(i,i+1,base4[b]);
        sfLFRbarcodeSnpMap[snp] = id;
      }
    }
  }
  //cout << "Finish loop in stlfr::ReadBarcodeList. sfLFRbarcodeMap size: " << sfLFRbarcodeMap.size() << endl;
  //cout << "sfLFRbarcodeSnpMap size: " << sfLFRbarcodeSnpMap.size() << endl;
}

void stlfr::findBarcode(Read* r) {
  int tried = 0, min = 0;
  int offsets[5] = {0, -1, 1, -2, 2};
  int stlfrPOS[3] = {mOptions->stlfr.pos1, mOptions->stlfr.pos2, mOptions->stlfr.pos3};
  int stlfrBarcodeLabels[3];
  //cout << "\nRun into stlfr::findBarcode. print stlfrPOS: " << stlfrPOS[0] << stlfrPOS[1] << stlfrPOS[2] << endl;
  //cout << "print sequence: " << r->mSeq.mStr << endl;

  //try to match traget sequences with known barcode list;
  while(tried < 5 && min == 0 ){
    int HitScore = 0;
    min = 2048;
    //int tryOffset = offsets[tried];
    for(int i = 0; i < 3; i++) {
      string thisCode = r->mSeq.mStr.substr(stlfrPOS[i] + offsets[tried], mOptions->stlfr.length);
      if(sfLFRbarcodeMap.count(thisCode) > 0){
        stlfrBarcodeLabels[i] = sfLFRbarcodeMap[thisCode];
        HitScore += 2* pow(10,i); //for debug
      }else if (sfLFRbarcodeSnpMap.count(thisCode) > 0){
        stlfrBarcodeLabels[i] = sfLFRbarcodeSnpMap[thisCode];
        HitScore += pow(10,i);
      }else {
        stlfrBarcodeLabels[i] = 0;
      }
      if (stlfrBarcodeLabels[i] < min)
        min = stlfrBarcodeLabels[i];
    }
    if(min == 0){
      missScore[offsets[tried]] += HitScore;
    }
    tried ++;
  }
  //for debug
  if(min > 0 ){
    tried --;
    //need to add an hash table for stats
  }else{
    //need to add an hash table for stats
  }
  //printf("%d\t%d\t%d\n",stlfrBarcodeLabels[0],stlfrBarcodeLabels[1],stlfrBarcodeLabels[2]);
  //sprintf(tag,"%8s",itoa(stlfrBarcodeLabels[0]));
  sprintf(stlfrBarcodeTag, "%04d_%04d_%04d",stlfrBarcodeLabels[0],stlfrBarcodeLabels[1],stlfrBarcodeLabels[2]);
  //printf("Print stlfrBarcodeTag: %s\n",stlfrBarcodeTag);
}

void stlfr::process(Read* r1, Read* r2) {
  if(!mOptions->stlfr.enabled)
      return;

  stringmap missScore, stlfrDebugs;
  int trys = 0, min = 0, hitScore = 0; //hitScore is used for debug stat

  if(mOptions->stlfr.loc == STLFR_LOC_READ1){
      findBarcode(r1);
      r1->keepFront(mOptions->stlfr.pos1);
  }
  else if(mOptions->stlfr.loc == STLFR_LOC_READ2 && r2){
      findBarcode(r2);
      r2->keepFront(mOptions->stlfr.pos1);
  }

  addstlfrToName(r1, stlfrBarcodeTag);
  if(r2){
      addstlfrToName(r2, stlfrBarcodeTag);
  }

}

void stlfr::addstlfrToName(Read* r, string stlfr){
  string tag;
  tag = "/" + stlfr;
  //string name = r->mName;
  //printf("Print  r->mName: %s\n",r->mName);
  int nameTailPos = r->mName.rfind("/1");
  if(nameTailPos == -1){
    nameTailPos = r->mName.rfind("/2");
  }
  if(nameTailPos == -1){
    nameTailPos = r->mName.find(' ');
  }
  r->mName = r->mName.replace(nameTailPos,0,tag);
}


bool stlfr::test() {
  return true;
}
