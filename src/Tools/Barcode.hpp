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
    ConstantBarcode(std::string inPattern, int inMismatches) : Barcode(inMismatches),pattern(inPattern) {}
    bool match_pattern(std::string sequence, int offset, int seq_start, int seq_end, fastqStats& stats);

    private:
    std::string pattern;
};
class VariableBarcode : public Barcode
{

    public:
    VariableBarcode(std::vector<std::string> inPatterns, int inMismatches) : Barcode(inMismatches),patterns(inPatterns) {}
    bool match_pattern(std::string sequence, int offset, int seq_start, int seq_end, fastqStats& stats);

    private:
    std::vector<std::string> patterns;

};