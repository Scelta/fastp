#include "processor.h"
#include "peprocessor.h"
#include "seprocessor.h"
#include "overlapanalysis.h"

Filter::Filter(Options* opt){
    mOptions = opt;
}


Filter::~Filter(){
}

int Filter::passFilter(Read* r) {
    if(r == NULL || r->length()==0) {
        return FAIL_LENGTH;
    }

    int rlen = r->length();
    int lowQualNum = 0;
    int nBaseNum = 0;

    // need to recalculate lowQualNum and nBaseNum if the corresponding filters are enabled
    if((mOptions->qualfilter.enabled && !mOptions->qualityCut.enabledOA) || mOptions->lengthFilter.enabled) {
        const char* seqstr = r->mSeq.mStr.c_str();
        const char* qualstr = r->mQuality.c_str();

        for(int i=0; i<rlen; i++) {
            char base = seqstr[i];
            char qual = qualstr[i];

            if(qual < mOptions->qualfilter.qualifiedQual)
                lowQualNum ++;

            if(base == 'N')
                nBaseNum++;
        }
    }

    if(mOptions->qualfilter.enabled && !mOptions->qualityCut.enabledOA) {
        if(lowQualNum > (mOptions->qualfilter.unqualifiedPercentLimit * rlen / 100.0) )
            return FAIL_QUALITY;
        else if(nBaseNum > mOptions->qualfilter.nBaseLimit )
            return FAIL_N_BASE;
    }

    if(mOptions->lengthFilter.enabled) {
        if(rlen < mOptions->lengthFilter.requiredLength)
            return FAIL_LENGTH;
        if(mOptions->lengthFilter.maxLength > 0 && rlen > mOptions->lengthFilter.maxLength)
            return FAIL_TOO_LONG;
    }

    if(mOptions->complexityFilter.enabled) {
        if(!passLowComplexityFilter(r))
            return FAIL_COMPLEXITY;
    }

    return PASS_FILTER;
}

bool Filter::passLowComplexityFilter(Read* r) {
    int diff = 0;
    int length = r->length();
    if(length <= 1)
        return false;
    const char* data = r->mSeq.mStr.c_str();
    for(int i=0; i<length-1; i++) {
        if(data[i] != data[i+1])
            diff++;
    }
    if( (double)diff/(double)(length-1) >= mOptions->complexityFilter.threshold )
        return true;
    else
        return false;
}

Read* Filter::trimAndCut(Read* r, int front, int tail) {
    // return the same read for speed if no change needed
    if(front == 0 && tail == 0 && !mOptions->qualityCut.enabled5 && !mOptions->qualityCut.enabled3 && !mOptions->qualityCut.enabledOA)
        return r;


    int rlen = r->length() - front - tail ;
    if (rlen < 0)
        return NULL;

    if(front == 0 && !mOptions->qualityCut.enabled5 && !mOptions->qualityCut.enabled3 && !mOptions->qualityCut.enabledOA){
        r->resize(rlen);
        return r;
    } else if(!mOptions->qualityCut.enabled5 && !mOptions->qualityCut.enabled3 && !mOptions->qualityCut.enabledOA){
        r->mSeq.mStr = r->mSeq.mStr.substr(front, rlen);
        r->mQuality = r->mQuality.substr(front, rlen);
        return r;
    }

    // need quality cutting

    int w = mOptions->qualityCut.windowSize;
    int l = r->length();
    const char* qualstr = r->mQuality.c_str();
    const char* seq = r->mSeq.mStr.c_str();
    // quality cutting forward
    if(mOptions->qualityCut.enabled5) {
        int s = front;
        if(l - front - tail - w <= 0)
            return NULL;

        int totalQual = 0;

        // preparing rolling
        for(int i=0; i<w-1; i++)
            totalQual += qualstr[s+i];

        for(s=front; s+w<l-tail; s++) {
            totalQual += qualstr[s+w-1];
            // rolling
            if(s > front) {
                totalQual -= qualstr[s-1];
            }
            // add 33 for phred33 transforming
            if((double)totalQual / (double)w >= 33 + mOptions->qualityCut.quality)
                break;
        }

        // the trimming in front is forwarded and rlen is recalculated
        if(s >0 )
            s = s+w-1;
        while(s<l && seq[s] == 'N')
            s++;
        front = s;
        rlen = l - front - tail;
    }

    // quality cutting backward
    if(mOptions->qualityCut.enabled3) {
        if(l - front - tail - w <= 0)
            return NULL;

        int totalQual = 0;
        int t = l - tail - 1;

        // preparing rolling
        for(int i=0; i<w-1; i++)
            totalQual += qualstr[t-i];

        for(t=l-tail-1; t-w>=front; t--) {
            totalQual += qualstr[t-w+1];
            // rolling
            if(t < l-tail-1) {
                totalQual -= qualstr[t+1];
            }
            // add 33 for phred33 transforming
            if((double)totalQual / (double)w >= 33 + mOptions->qualityCut.quality)
                break;
        }

        if(t < l-1)
            t = t-w+1;
        while(t>=0 && seq[t] == 'N')
            t--;
        rlen = t - front + 1;
    }

    // quality cutting by Overall Accuracy
    if(mOptions->qualityCut.enabledOA) {
        int qualScore[rlen];
        //Table obtained by lg(1-10^(-PhedQuality/10))*1e10)
        const double Accstr[50] = {
          -1.0000000, -0.6868253, -0.4329234, -0.3020624, -0.2204808,
          -0.1650885, -0.1256276, -0.0966529, -0.0749404, -0.0584352,
          -0.0457575, -0.0359445, -0.0283048, -0.0223307, -0.0176431,
          -0.0139554, -0.0110483, -0.0087529, -0.0069382, -0.0055022,
          -0.0043648, -0.0034635, -0.0027489, -0.0021821, -0.0017324,
          -0.0013755, -0.0010923, -0.0008674, -0.0006889, -0.0005471,
          -0.0004345, -0.0003451, -0.0002741, -0.0002177, -0.0001729,
          -0.0001374, -0.0001091, -0.0000867, -0.0000688, -0.0000547,
          -0.0000434, -0.0000345, -0.0000274, -0.0000218, -0.0000173,
          -0.0000137, -0.0000109, -0.0000087, -0.0000069, -0.0000055};

        for(int i=0;i<rlen;i++)
            qualScore[i]=qualstr[i]-33;

        double OAseedAcc = Accstr[mOptions->qualityCut.OAsqual];
        double OAfragAcc = Accstr[mOptions->qualityCut.OAfqual];
        int s = front;
        int the_s = s;
        if(l - front - tail - w <= 0)
            return NULL;

        double seedAcc[l - front - tail - w];
        seedAcc[0] = 0;

        // preparing 1st seed
        for(int i=0; i<w; i++)
            seedAcc[s] += Accstr[qualScore[s+i]];

        // Rolling to find a better seed

        while(seedAcc[s] < OAseedAcc && s+w<l-tail){
            s++;
            seedAcc[s] = seedAcc[s-1] + Accstr[qualScore[s+w-1]] - Accstr[qualScore[s-1]];
            the_s = (seedAcc[s] > seedAcc[the_s])?s:the_s;
        }

        // Extend the fragment's tail
        double fragAcc = seedAcc[the_s];
        int f=the_s;
        double ignore_the_lowest_acc = 0;
        double theFragAcc = fragAcc;

        while(f+w<l-tail){
            if(Accstr[qualScore[f+w]] < ignore_the_lowest_acc){
                theFragAcc += ignore_the_lowest_acc;
                ignore_the_lowest_acc = Accstr[qualScore[f+w]];
            }else{
                theFragAcc += Accstr[qualScore[f+w]];
            }

            if (theFragAcc <= OAfragAcc)
                break;

            f++;
            fragAcc = theFragAcc;
        }
        // Discard if accuracy is too low
        if(f+1+w >= l-tail && fragAcc < OAfragAcc)
            return NULL;
        // Dtermine the front and length
        front = the_s;
        rlen = f - the_s + w + 1;
    }

    if(rlen <= 0 || front >= l-1)
        return NULL;

    r->mSeq.mStr = r->mSeq.mStr.substr(front, rlen);
    r->mQuality = r->mQuality.substr(front, rlen);

    return r;
}

bool Filter::filterByIndex(Read* r) {
    if(mOptions->indexFilter.enabled) {
        if( match(mOptions->indexFilter.blacklist1, r->firstIndex(), mOptions->indexFilter.threshold) )
            return true;
    }
    return false;
}

bool Filter::filterByIndex(Read* r1, Read* r2) {
    if(mOptions->indexFilter.enabled) {
        if( match(mOptions->indexFilter.blacklist1, r1->firstIndex(), mOptions->indexFilter.threshold) )
            return true;
        if( match(mOptions->indexFilter.blacklist2, r2->lastIndex(), mOptions->indexFilter.threshold) )
            return true;
    }
    return false;
}

bool Filter::match(vector<string>& list, string target, int threshold) {
    for(int i=0; i<list.size(); i++) {
        int diff = 0;
        int len1 = list[i].length();
        int len2 = target.length();
        for(int s=0; s<len1 && s<len2; s++) {
            if(list[i][s] != target[s]) {
                diff++;
                if(diff>threshold)
                    break;
            }
        }
        if(diff <= threshold)
            return true;
    }
    return false;
}

bool Filter::test() {
    Read r("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTT",
        "+",
        "/////CCCCCCCCCCCC////CCCCCCCCCCCCCC////E");
    Options opt;
    opt.qualityCut.enabled5 = true;
    opt.qualityCut.enabled3 = true;
    opt.qualityCut.windowSize = 4;
    opt.qualityCut.quality = 20;
    Filter filter(&opt);
    Read* ret = filter.trimAndCut(&r, 0, 1);
    ret->print();

    return ret->mSeq.mStr == "CCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        && ret->mQuality == "CCCCCCCCCCC////CCCCCCCCCCCCC";
}

bool Filter::testOA() {
    Read r("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTT",
        "+",
        "/////CCCCCCCCCCCC////CCCCCCCCCCCCCC////E");
    Options opt;
    opt.qualityCut.enabled5 = false;
    opt.qualityCut.enabled3 = false;
    opt.qualityCut.enabledOA = true;
    opt.qualityCut.windowSize = 10;
    opt.qualityCut.OAsqual = 20;
    opt.qualityCut.OAfqual = 10;
    Filter filter(&opt);
    Read* ret = filter.trimAndCut(&r, 0, 1);
    ret->print();

    return ret->mSeq.mStr == "ACCCCCCCCCCCCCCC"
        && ret->mQuality  == "CCCCCCCCCCCC////";
}
