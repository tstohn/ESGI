#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <thread>

//new datatypes
class Barcode
{
    public:
    Barcode(int inMismatches) : mismatches(inMismatches) {}

    private:
    int mismatches;

    //overwritten function to match sequence pattern(s)
    virtual std::string match_pattern(std::string sequence) = 0;

};

class ConstantBarcode : public Barcode
{

    public:
    ConstantBarcode(std::string inSequence, int inMismatches) : Barcode(inMismatches),pattern(inSequence) {}
    std::string match_pattern();

    private:
    std::string pattern;
};
class VariableBarcode : public Barcode
{

    public:
    VariableBarcode(std::vector<std::string> inSequences, int inMismatches) : Barcode(inMismatches),patterns(inSequences) {}
    std::string match_pattern();

    private:
    std::vector<std::string> patterns;

};

typedef std::shared_ptr<Barcode> barcodePtr;
typedef std::vector<barcodePtr> barcodePtrVector; 