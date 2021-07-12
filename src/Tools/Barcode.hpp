#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <thread>

#include "helper.hpp"

class Barcode;
typedef std::shared_ptr<Barcode> BarcodePatternPtr;
typedef std::vector<BarcodePatternPtr> BarcodePatternVector; 
typedef std::shared_ptr<BarcodePatternVector> BarcodePatternVectorPtr; 

//new datatypes
class Barcode
{
    public:
    Barcode(int inMismatches) : mismatches(inMismatches) {}
    int mismatches;

    private:
    //overwritten function to match sequence pattern(s)
    virtual bool match_pattern(std::string sequence, int offset, int seq_start, int seq_end, fastqStats& stats) = 0;

};

class ConstantBarcode : public Barcode
{

    public:
    ConstantBarcode(std::string inPattern, int inMismatches) : pattern(inPattern),Barcode(inMismatches) {}
    bool match_pattern(std::string sequence, int offset, int seq_start, int seq_end, fastqStats& stats)
    {
    int start = 0, end = 0, score = 0;

    sequence.erase(0, offset);
    // start in seq is at: start-1, end-start+1
    if(levenshtein(sequence, pattern, mismatches, start, end, score))
    {
        int startIdx = start-1;
        int endIdx = end; // inclusive
        //minor mismatches that are allowed per mb and ma
        ++stats.moderateMatches;
        return true;
    }
    else
    {
        //bigger mismatches
        ++stats.noMatches;
        return false;
    }

}

    private:
    std::string pattern;
};
class VariableBarcode : public Barcode
{

    public:
    VariableBarcode(std::vector<std::string> inPatterns, int inMismatches) : patterns(inPatterns), Barcode(inMismatches) {}
    bool match_pattern(std::string sequence, int offset, int seq_start, int seq_end, fastqStats& stats)
    {
        return true;
    }

    private:
    std::vector<std::string> patterns;

};