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
        const double Accstr[50] = {
          -1.00000000, -0.68682512, -0.43292312, -0.30206212, -0.22048112, -0.16508912, -0.12562812, -0.09665312, -0.07494012, -0.05843512,
          -0.04575712, -0.03594512, -0.02830512, -0.02233112, -0.01764312, -0.01395512, -0.01104812, -0.00875312, -0.00693812, -0.00550212,
          -0.00436512, -0.00346312, -0.00274912, -0.00218212, -0.00173212, -0.00137612, -0.00109212, -0.00086712, -0.00068912, -0.00054712,
          -0.00043512, -0.00034512, -0.00027412, -0.00021812, -0.00017312, -0.00013712, -0.00010912, -0.00008712, -0.00006912, -0.00005512,
          -0.00004312, -0.00003412, -0.00002712, -0.00002212, -0.00001712, -0.00001412, -0.00001112, -0.00000912, -0.00000712, -0.00000512};

        double OAseedAcc = Accstr[mOptions->qualityCut.OAsqual];
        double OAfragAcc = Accstr[mOptions->qualityCut.OAfqual];
        int s = front;
        int the_s = s;
        if(l - front - tail - w <= 0)
            return NULL;

        double seedAcc[l - front - tail - w];

        // preparing 1st seed
        for(int i=0; i<w-1; i++)
            seedAcc[s] += Accstr[qualstr[s+i]];

        // Rolling to find a better seed
        while(seedAcc[s] < OAseedAcc && s+w<l-tail){
            s++;
            seedAcc[s] = seedAcc[s-1] + Accstr[qualstr[s+w-1]] - Accstr[qualstr[s-1]];
            the_s = (seedAcc[s] > seedAcc[s-1])?s:s-1;
        }

        // Extend the fragment's tail
        double fragAcc = seedAcc[the_s];
        int f=the_s;
        double ignore_the_lowest_acc = 0;
        while(fragAcc > OAfragAcc && f+w<l-tail){
            f++;
            if(Accstr[qualstr[f+w-1]] < ignore_the_lowest_acc){
                fragAcc += ignore_the_lowest_acc;
                ignore_the_lowest_acc = Accstr[qualstr[f+w-1]];
            }else{
                fragAcc += Accstr[qualstr[f+w-1]];
            }
        }
        // Discard if accuracy is too low
        if(fragAcc < OAfragAcc)
            return NULL;
        // Dtermine the front and length
        front = the_s;
        rlen = f - s + w;
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
