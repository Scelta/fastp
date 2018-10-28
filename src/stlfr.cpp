#include "stlfr.h"

stringmap sfLFRbarcodeMap, sfLFRbarcodeSnpMap;
intmap missScore, stlfrLog;

stlfr::stlfr(Options* opt){
    mOptions = opt;
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
  mOptions->stlfr.barcodeSpace = sfLFRbarcodeMap.size();
  #define BSIZE sfLFRbarcodeMap.size();
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
  r->stlfrB1 = stlfrBarcodeLabels[0];
  r->stlfrB2 = stlfrBarcodeLabels[1];
  r->stlfrB3 = stlfrBarcodeLabels[2];
  char stlfrBarcodeTag[16];
  sprintf(stlfrBarcodeTag, "%04d_%04d_%04d",stlfrBarcodeLabels[0],stlfrBarcodeLabels[1],stlfrBarcodeLabels[2]);
  addstlfrToName(r, stlfrBarcodeTag);
  printf("Print stlfrBarcodeTag: %s\n",stlfrBarcodeTag);
}

void stlfr::process(Read* r1, Read* r2) {
  if(!mOptions->stlfr.enabled)
      return;

  stringmap missScore, stlfrDebugs;
  int trys = 0, min = 0, hitScore = 0; //hitScore is used for debug stat
  string thisBarcodes;
  if(mOptions->stlfr.loc == STLFR_LOC_READ1){
      findBarcode(r1);
      thisBarcodes = r1->mSTLFR;
      r1->keepFront(mOptions->stlfr.pos1);
  }
  else if(mOptions->stlfr.loc == STLFR_LOC_READ2 && r2){
      findBarcode(r2);
      thisBarcodes = r2->mSTLFR;
      addstlfrToName(r1, thisBarcodes);
      r2->keepFront(mOptions->stlfr.pos1);
  }
}

void stlfr::addstlfrToName(Read* r, string stlfr){
  string tag;
  tag = "/" + stlfr;
  r->mSTLFR = stlfr;
  int nameTailPos = r->mName.rfind("/1");
  if(nameTailPos == -1){
    nameTailPos = r->mName.rfind("/2");
  }
  if(nameTailPos == -1){
    nameTailPos = r->mName.find(' ');
  }
  r->mName = r->mName.replace(nameTailPos,0,tag);
  //cout << "add tagName to : " << r->mName << endl;

}


bool stlfr::test() {
  return true;
}


//init for stlfr stats
void stlfr::stlfrStats(Options* opt){
  int sbSize = mOptions->stlfr.barcodeSpace + 1;
  mStlfrBarcode = new int **[sbSize];
  for(int i=0; i<=sbSize; i++){
    mStlfrBarcode[i] = new int *[sbSize];
    for(int j=0; j<=sbSize; i++){
      mStlfrBarcode[i][j] = new int [sbSize];
    }
  }
  memset(mStlfrBarcode, 0, sizeof(long)*sbSize*sbSize*sbSize);
  cout << "init stat" << endl;
  mStlfrValid = 0;
}

stlfr::~stlfr(){
  //delete array
  int sbSize = mOptions->stlfr.barcodeSpace + 1;
  for(int i=0; i<=sbSize; i++){
    for(int j=0; j<=sbSize; i++){
      delete mStlfrBarcode[i][j];
    }
    delete mStlfrBarcode[i];
  }
  delete mStlfrBarcode;
}

void stlfr::statStlfr(Read* r){
  int b1 = r->stlfrB1;
  int b2 = r->stlfrB2;
  int b3 = r->stlfrB3;
  mStlfrBarcode[b1][b2][b3] ++;
  cout << "check" << endl;
  printf("stat stlfr: %d_%d_%d:%d\n",b1,b2,b3,mStlfrBarcode[b1][b2][b3]);
}

//merge stlfr barcode frequency
stlfr* stlfr::merge(vector<stlfr*>& list) {
  if(list.size() == 0)
      return NULL;

  // barcode stats
  stlfr* s = new stlfr(list[0]->mOptions);
  int sRange = BSIZE;
  for(int t=0; t<list.size(); t++) {
    for(int b1=0; b1<=sRange; b1++) {
      for(int b2=0; b2<=sRange; b2++) {
        for(int b3=0; b3<=sRange; b3++) {
          if(list[t]->mStlfrBarcode[b1][b2][b3] >0 && s->mStlfrBarcode[b1][b2][b3]==0 && b1 != 0 && b2 != 0 && b3 != 0)
            s->mStlfrValid ++;
          s->mStlfrBarcode[b1][b2][b3] += list[t]->mStlfrBarcode[b1][b2][b3];
        }
      }
    }
  }
  return s;
}

void stlfr::print() {
    cerr << "Valid stlfr barcodes found: " << mStlfrValid << endl;
}

void stlfr::reportJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    // barcode frequency
    ofs << padding << "\t" << "\"stLFR_barcode_frequency\": {" << endl;
    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    for(int b1=0; b1<=barcodeSpace; b1++) {
      count1 = 0;
      for(int b2=0; b2<=barcodeSpace; b2++) {
        count2 = 0;
        for(int b3=0; b3<=barcodeSpace; b3++) {
          count3 = 0;
          long count = mStlfrBarcode[b1][b2][b3];
          if(mStlfrBarcode[b1][b2][b3] != 0){
            ofs << padding << "\t\"" << b1 << "_" << b2 << "_" << b3 << "\":" << count;
            count3 ++;
          }
        }
        count2 += count3;
      }
      if(count2 !=0){
        ofs << "," << endl;
      }
    }
    ofs << padding << "\t" << "}," << endl;

    ofs << padding << "}," << endl;
}
