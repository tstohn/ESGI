#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <thread>
#include <unordered_set>

#include "helper.hpp"

class Barcode;
typedef std::shared_ptr<Barcode> BarcodePtr;
typedef std::vector<BarcodePtr> BarcodeVector; 
typedef std::shared_ptr<BarcodeVector> BarcodeVectorPtr; 

//class to handle a barcode pattern 
//can be used to iterate through the barcode, and stores additional information like:
//it stores if the pattern contains DNA barcodes which require different handling
class BarcodePattern
{
    public:
        //constructor
        BarcodePattern(bool dna, std::string name, BarcodeVectorPtr pattern) :  
        containsDNA(dna), patternName(name), barcodePattern(pattern) {}

        //class variables
        bool containsDNA;
        std::string patternName; //this is also the file this pattern will be written to
        BarcodeVectorPtr barcodePattern;

        //class functions
        //write multiplexed lines to file (must be specific for DNA, AB-barcodes, etc.)
        void write_demultiplexed_line(const std::vector<std::string> barcodeList, std::string dna = "");

        //barcodeVector/ iterator functions
        // Add a barcode to the barcodePattern
        void add_barcode(const BarcodePtr& barcode) {
            barcodePattern->push_back(barcode);
        }
        // Get the size of the barcodePattern
        std::size_t size() const {
            return barcodePattern->size();
        }
        // Access element by index
        BarcodePtr& operator[](std::size_t index) {
            return (*barcodePattern)[index];
        }
        const BarcodePtr& operator[](std::size_t index) const {
            return (*barcodePattern)[index];
        }
        // Iterator types
        using iterator = typename std::vector<BarcodePtr>::iterator;
        using const_iterator = typename std::vector<BarcodePtr>::const_iterator;
        using reverse_iterator = typename std::vector<BarcodePtr>::reverse_iterator;
        using const_reverse_iterator = typename std::vector<BarcodePtr>::const_reverse_iterator;
        // Begin and end iterators
        iterator begin() {
            return barcodePattern->begin();
        }
        const_iterator begin() const {
            return barcodePattern->begin();
        }
        iterator end() {
            return barcodePattern->end();
        }
        const_iterator end() const {
            return barcodePattern->end();
        }
        // Reverse iterators
        reverse_iterator rbegin() {
            return barcodePattern->rbegin();
        }
        const_reverse_iterator rbegin() const {
            return barcodePattern->rbegin();
        }
        reverse_iterator rend() {
            return barcodePattern->rend();
        }
        const_reverse_iterator rend() const {
            return barcodePattern->rend();
        }
};

typedef std::shared_ptr<BarcodePattern> BarcodePatternPtr; 
typedef std::shared_ptr<std::vector<BarcodePatternPtr>> MultipleBarcodePatternVectorPtr; 

struct mappingSolution{             
    int seq_start;
    int seq_end;
    int score = INT_MAX;
    std::string realBarcode;
    int differenceInBarcodeLength;
}; 

//new datatypes
class Barcode
{
    public:
    //per default the length of a barcode is set to 0 (unknown), only for UMIs it must be set
    Barcode(std::string name, int inMismatches, int inLength = 0) : name(name), mismatches(inMismatches), length(inLength) {}
    //virtual destructor, needed to be called when destructing classes that inherit from Barcode
    virtual ~Barcode() = default;

    std::string name;
    int mismatches;
    unsigned int length; //length of zero means the length of this barcode is unknown

    //reverse complement is a Barcode function that should be available globally
    static std::string generate_reverse_complement(std::string seq)
    {
        auto lambda = [](const char c) 
        {
            switch (c) 
            {
                case 'A':
                    return 'T';
                case 'G':
                    return 'C';
                case 'C':
                    return 'G';
                case 'T':
                    return 'A';
                case 'N':
                    return 'N';
                default:
                    throw std::domain_error("Invalid nucleotide.");
            }
        };
        std::string newSeq = seq;
        std::transform(seq.crbegin(), seq.crend(), newSeq.begin(), lambda);
        return newSeq;
    }
    //overwritten function to match sequence pattern(s)
    virtual bool align(std::string& matchedBarcode, const std::string& fastqLine,const int targetOffset,
                       int& targetEnd, int& delNum, int& insNum, int& substNum,
                       bool reverse = false) = 0;
    virtual std::vector<std::string> get_patterns() = 0;
    virtual bool is_wildcard() = 0;
    virtual bool is_constant() = 0;
    virtual bool is_stop() = 0;
    virtual bool is_dna() = 0;
    virtual bool is_read_end() = 0;

};

class ConstantBarcode : public Barcode
{

    public:
    ConstantBarcode(std::string inPattern, int inMismatches) : Barcode(inPattern, inMismatches), pattern(inPattern)
    {
        revCompPattern = generate_reverse_complement(pattern);
        // Configure Edlib
        config = edlibNewAlignConfig(
            inMismatches,        // Maximum allowed edit distance
            EDLIB_MODE_SHW,     // Semi-global alignment (deletions in the target at the end are not penalized), 
                                // the developers refer to it also as 'Prefix method'
            EDLIB_TASK_PATH,    // Request full alignment path (M, I, D, S)
            NULL, 0             // No custom alphabet
        );
    }

    bool align(std::string& matchedBarcode, const std::string& fastqLine,const int targetOffset,
               int& targetEnd, int& delNum, int& insNum, int& substNum,
               bool reverse = false)
    {
        bool foundAlignment = false;

        //get length of substring, length does not depend on reverse/ forward pattern
        std::string target;
        int substringLength = pattern.length()+mismatches;
        if(targetOffset + substringLength > fastqLine.size()){substringLength = fastqLine.size()-targetOffset;};
        //std::cout << "\t LENGTH: " <<substringLength << " seq: " << fastqLine.size()<< "\n";
        target = fastqLine.substr(targetOffset, substringLength);

        //set the pattern to use for reverse or forward mapping
        std::string usedPattern = pattern;
        if(reverse){usedPattern = revCompPattern;}

        //map the pattern to the target sequence
        foundAlignment = run_alignment(usedPattern, target, targetEnd, config, delNum,  insNum, substNum);
        matchedBarcode = pattern;

        return foundAlignment;
    }

    std::vector<std::string> get_patterns()
    {
        std::vector<std::string> patterns = {pattern};
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return true;}
    bool is_stop(){return false;}
    bool is_dna(){return false;}
    bool is_read_end(){return false;}

    private: 
        std::string pattern;
        std::string revCompPattern;
        EdlibAlignConfig config;
};
class VariableBarcode : public Barcode
{

    public:
    VariableBarcode(std::vector<std::string> inPatterns, std::string name, int inMismatches) : Barcode(name, inMismatches), patterns(inPatterns)
    {
        for(std::string pattern : patterns)
        {
            std::string revCompPattern = generate_reverse_complement(pattern);
            revCompPatterns.push_back(revCompPattern);
        }

        config = edlibNewAlignConfig(
            inMismatches,        // Maximum allowed edit distance
            EDLIB_MODE_SHW,     // Semi-global alignment (deletions in the target at the end are not penalized), 
                                // the developers refer to it also as 'Prefix method'
            EDLIB_TASK_PATH,    // Request full alignment path (M, I, D, S)
            NULL, 0             // No custom alphabet
        );

        //instant look-up table for perfect barcodes
        equalLengthBarcodes = true;
        size_t lengthOne = patterns.at(0).size();
        for (const std::string& pattern : patterns) 
        {
            if (pattern.size() != lengthOne) 
            {
                equalLengthBarcodes = false;
                break;
            }
        }
        if(equalLengthBarcodes)
        {
            barcodeSet.insert(patterns.begin(), patterns.end());
        }

        //calculate minimum conversion rates of barcodes

    }

    bool align(std::string& matchedBarcode, const std::string& fastqLine, const int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)
    {
        //check if we can instantly match pattern
        if(equalLengthBarcodes && (fastqLine.size() >= (targetOffset + patterns.at(0).size())))
        {
            std::string target = fastqLine.substr(targetOffset, patterns.at(0).size());
            if (barcodeSet.find(target) != barcodeSet.end()) 
            {
                targetEnd = patterns.at(0).size();
                matchedBarcode = target;
                return true;
            }  
        }

        std::vector<std::string> patternsToMap = patterns;
        //get the reverse pattern list if we have reverse string
        if(reverse){patternsToMap = revCompPatterns;}

        std::string bestFoundPattern;
        bool bestFoundAlignment = false;
        int bestTargetEnd = -1;
        int bestEditDist = mismatches+1;
        for(int patternIdx = 0; patternIdx!= patternsToMap.size(); ++patternIdx)
        {
            bool foundAlignment = false;

            int delNumTmp;
            int insNumTmp;
            int substNumTmp;
            delNumTmp=insNumTmp=substNumTmp=targetEnd=0;
            //get pattern, its length can vary
            std::string usedPattern = patternsToMap.at(patternIdx);

            //define target sequence (can differ for every barcode due to its length)
            std::string target;
            int substringLength = usedPattern.length()+mismatches;
            //std::cout << "\t LENGTH: " <<substringLength << " seq: " << fastqLine.size()<< "\n";
            if(targetOffset + substringLength > fastqLine.size()){substringLength = fastqLine.size()-targetOffset;};
            target = fastqLine.substr(targetOffset, substringLength);

            //std::cout << " in barcode trying barcode: " << usedPattern << " with target sequ " << target<< "\n";

            //map the pattern to the target sequence
            foundAlignment = run_alignment(usedPattern, target, targetEnd, config, delNumTmp,  insNumTmp, substNumTmp);
            
            if(foundAlignment && (delNum+insNum+substNum)<bestEditDist)
            {
                bestFoundAlignment = foundAlignment;
                bestFoundPattern = patterns.at(patternIdx); //the barcode is the TRUE forward barcode, not the reverse complement
                bestTargetEnd = targetEnd;
                bestEditDist = (delNum+insNum+substNum);
                delNum = delNumTmp;
                insNum = insNumTmp;
                substNum = substNumTmp;
            }

            if(bestEditDist<=1)
            {
                break;
            }
        }

        targetEnd = bestTargetEnd;
        matchedBarcode = bestFoundPattern;

        return bestFoundAlignment;
    }

    std::vector<std::string> get_patterns()
    {
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return false;}
    bool is_stop(){return false;}
    bool is_dna(){return false;}
    bool is_read_end(){return false;}

    private:
        std::vector<std::string> patterns;
        std::vector<std::string> revCompPatterns;
        std::unordered_set<std::string> barcodeSet;
        bool equalLengthBarcodes;

        EdlibAlignConfig config;
};

class WildcardBarcode : public Barcode
{
    // TODO: move public parameter mismatches to a private and derived parameter
    //wildcardBarcode doe snot make use of mismatches yet, since anyways we do not know the sequence,
    //therefore its an unused parameter, just set for completeness as these classes derive from Barcode (initialized with mismatches, see up...)
    public:
    WildcardBarcode(int inMismatches, std::string name, int inLength) : Barcode(name, inMismatches, inLength) {}

    bool align(std::string& matchedBarcode, const std::string& target, const int positionInFastqLine,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)
    {
        matchedBarcode = target.substr(positionInFastqLine, length);
        targetEnd = (target.length() < length) ? target.length() : length;
        // e.g.: [AGTAGT]cccc: start=0 end=6 end is first not included idx
        return true;
    }

    std::vector<std::string> get_patterns()
    {
        std::vector<std::string> patterns = {std::string(length, 'X')};
        return patterns;
    }
    bool is_wildcard(){return true;}
    bool is_constant(){return false;}
    bool is_stop(){return false;}
    bool is_dna(){return false;}
    bool is_read_end(){return false;}

};

//A stop-barcode: mapping on both sides is only done up to here
//the name ('*') is not writen into the output file
class StopBarcode : public Barcode
{
    public:
    StopBarcode(std::string inPattern, int inMismatches) : Barcode("*", inMismatches),pattern(inPattern) {}

    bool align(std::string& matchedBarcode, const std::string& target, const int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)

        {
            return false;
        }
    std::vector<std::string> get_patterns()
    {
        std::vector<std::string> patterns = {pattern};
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return false;}
    bool is_stop(){return true;}
    bool is_dna(){return false;}
    bool is_read_end(){return false;}

    private:
    std::string pattern; //just a string of "XXXXX"
};

//A stop-barcode: mapping on both sides is only done up to here
//the name ('*') is not writen into the output file
class ReadSeperatorBarcode : public Barcode
{
    public:
    ReadSeperatorBarcode(std::string inPattern, int inMismatches) : Barcode("-", inMismatches),pattern(inPattern) {}

    bool align(std::string& matchedBarcode, const std::string& target, const int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)

        {
            return false;
        }
    std::vector<std::string> get_patterns()
    {
        std::vector<std::string> patterns = {pattern};
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return false;}
    bool is_stop(){return false;}
    bool is_dna(){return false;}
    bool is_read_end(){return true;}

    private:
    std::string pattern; //just a string of "XXXXX"
};

//barcode representing a genomic region. We can set the number of mismatches or when set to -1 simply run it with default
//star settings
//pattern is simply "DNA"
class DNABarcode : public Barcode
{
    public:
    DNABarcode(int inMismatches = -1) : Barcode("DNA", inMismatches) {pattern = "DNA";}

    bool align(std::string& matchedBarcode, const std::string& target, const int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)

        {
            return false;
        }
    std::vector<std::string> get_patterns()
    {
        std::vector<std::string> patterns = {pattern};
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return false;}
    bool is_stop(){return false;}
    bool is_dna(){return true;}
    bool is_read_end(){return false;}

    private:
    std::string pattern; //just a string of "DNA"
};