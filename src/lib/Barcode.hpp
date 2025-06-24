#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <thread>
#include <unordered_set>
#include <unordered_map>

#include "helper.hpp"

class Barcode;
typedef std::shared_ptr<Barcode> BarcodePtr;
typedef std::vector<BarcodePtr> BarcodeVector; 
typedef std::shared_ptr<BarcodeVector> BarcodeVectorPtr; 

enum class PatternType 
{
    Forward,
    Reverse
};
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
        BarcodeVectorPtr detachedReversePattern; //in case we do not have one long pattern with a fw&rv read
        //but more two independent read that should be mapped seperately

        //class functions
        //write multiplexed lines to file (must be specific for DNA, AB-barcodes, etc.)
        void write_demultiplexed_line(const std::vector<std::string> barcodeList, std::string dna = "");

        //barcodeVector/ iterator functions
        // Add a barcode to the barcodePattern
        void add_barcode(const BarcodePtr& barcode, PatternType type = PatternType::Forward) 
        {
            get_pattern(type)->push_back(barcode);
        }
        // Get the size of the barcodePattern
        std::size_t size(PatternType type = PatternType::Forward) const 
        {
            return get_pattern(type)->size();
        }

        // Iterator types
        using iterator = typename std::vector<BarcodePtr>::iterator;
        using const_iterator = typename std::vector<BarcodePtr>::const_iterator;
        using reverse_iterator = typename std::vector<BarcodePtr>::reverse_iterator;
        using const_reverse_iterator = typename std::vector<BarcodePtr>::const_reverse_iterator;
        // Begin and end iterators
        iterator begin(PatternType type = PatternType::Forward) {
            return get_pattern(type)->begin();
        }
        const_iterator begin(PatternType type = PatternType::Forward) const {
            return get_pattern(type)->begin();
        }
        iterator end(PatternType type = PatternType::Forward) {
            return get_pattern(type)->end();
        }
        const_iterator end(PatternType type = PatternType::Forward) const {
            return get_pattern(type)->end();
        }
        // Reverse iterators
        reverse_iterator rbegin(PatternType type = PatternType::Forward) {
            return get_pattern(type)->rbegin();
        }
        const_reverse_iterator rbegin(PatternType type = PatternType::Forward) const {
            return get_pattern(type)->rbegin();
        }
        reverse_iterator rend(PatternType type = PatternType::Forward) {
            return get_pattern(type)->rend();
        }
        const_reverse_iterator rend(PatternType type = PatternType::Forward) const {
            return get_pattern(type)->rend();
        }

    private:
        BarcodeVectorPtr get_pattern(PatternType type) const 
        {
            if(type == PatternType::Reverse){return detachedReversePattern;}
            return barcodePattern; //if not reverse return forward pattern
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
    virtual bool align(std::string& matchedBarcode, const std::string& fastqLine, const unsigned int targetOffset,
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

    bool align(std::string& matchedBarcode, const std::string& fastqLine,const unsigned int targetOffset,
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
        calculate_barcode_conversionRates();

    }

    void calculate_barcode_conversionRates()
    {
        EdlibAlignConfig barcodeConversionConfig = edlibNewAlignConfig(
            -1,                 // no limit for edit distancve
            EDLIB_MODE_SHW,     // Semi-global alignment (deletions in the target at the end are not penalized), 
                                // the developers refer to it also as 'Prefix method', e.g. one barcode is ATC and another ATCATC
                                //this distance should still be zero...
            EDLIB_TASK_PATH,    // Request full alignment path (M, I, D, S)
            NULL, 0);             // No custom alphabet

        for (size_t i = 0; i < patterns.size(); ++i) 
        {
            const std::string a = patterns.at(i);
            int min_rate = std::numeric_limits<int>::max();
            int minElement = std::numeric_limits<int>::max();

            for (size_t j = 0; j < patterns.size(); ++j) 
            {
                if(i == j){continue;}

                const std::string b = patterns.at(j);
                int del = 0;
                int ins = 0;
                int subst = 0;
                int targetEnd;

                //for stagger barcode the first barcode has to be pattern, the secone the target, bcs. we do not cound deletions on the target
                // e.g.: barcode A and AGT should have a conversion rate of 0! So for staggered barcodes we need to test ALL barcodes and then take
                //the best match with the longest matching sequence...
                if (a.length() > b.length()) 
                {
                    run_alignment(b, a, targetEnd, barcodeConversionConfig, del, ins, subst);
                } else 
                {
                    run_alignment(a, b, targetEnd, barcodeConversionConfig, del, ins, subst);
                }

                int rate = del + ins + subst;
                if (rate < min_rate) 
                {
                    min_rate = rate;
                    minElement = j;
                }
            }

            //warning if barcodes are the same/ or one is suffix of the other
            if(min_rate == 0)
            {
                std::cout << "WARNING: The data contains barcodes with a mismatch distance of 0!!! Barcodes: " << patterns.at(i) << ", "<< patterns.at(minElement) << "\n" <<
                "In case these are barcodes of variable length, we map the longest barcode mapping with no errors!!\n";
            }

            //warning if minimum conversion is lower than allowed mismatches
            if(min_rate <= mismatches)
            {
                std::cout << "WARNING: \nFor barcodes in " << name << ": " << mismatches << " mismatches are allowed, but with " << min_rate <<
                " mismatches we can already convert " << patterns.at(i)  << " into "<< patterns.at(minElement)  <<  " (semi-global alignment with un-punished deletions in target)!!! We can still find a best-fitting barcode, but you might reconsider the choice of allowed mistmaches!\n";
            }

            pattern_conversionrates[a] = min_rate;

        }
    }

    bool align(std::string& matchedBarcode, const std::string& fastqLine, const unsigned int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)
    {
        //check if we can instantly match pattern
        if(equalLengthBarcodes && (fastqLine.size() >= (targetOffset + patterns.at(0).size())))
        {
            std::string target = fastqLine.substr(targetOffset, patterns.at(0).size());
            if(reverse){target = generate_reverse_complement(target);}
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
        bool severalMatches = false; //if there are several best-fitting solutions (only the case when the number of allowed mismatches
        //is bigger than possible barcode-conversion numbers) we discard the solution

        for(size_t patternIdx = 0; patternIdx!= patternsToMap.size(); ++patternIdx)
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

            //map the pattern to the target sequence
            foundAlignment = run_alignment(usedPattern, target, targetEnd, config, delNumTmp,  insNumTmp, substNumTmp);
            
            if(foundAlignment && (delNumTmp+insNumTmp+substNumTmp)<=bestEditDist)
            {
                //if we have variable length barcodes, and two barcodes map equally well, 
                //we store the barcode of a longer sequence
                //example: barcodes = ["ATC", "ATCATC"], if we find the barcode ATCATC we store that one and not ATC...
                //this if clause checks if the new match is shorter than odl one, if so do not store it
                if(!equalLengthBarcodes && targetEnd < bestTargetEnd)
                {
                    continue;
                }

                //check if we have several best solutions
                if((delNumTmp+insNumTmp+substNumTmp)==bestEditDist)
                {
                    severalMatches = true;
                }
                else
                {
                    severalMatches = false;
                }

                bestFoundAlignment = foundAlignment;
                bestFoundPattern = patterns.at(patternIdx); //the barcode is the TRUE forward barcode, not the reverse complement
                bestTargetEnd = targetEnd;
                bestEditDist = (delNumTmp+insNumTmp+substNumTmp);
                delNum = delNumTmp;
                insNum = insNumTmp;
                substNum = substNumTmp;

                //if we found a new best match, check if this is already the best match we can ever get (minimal conversion dist between barcodes)
                int minConversion = pattern_conversionrates.at(patterns.at(patternIdx));
                //we can do this since levenshtein distance fullfills the triangle inequality is a distance metric
                //imagine there is a second barcode that could fit better: this second barcode must have a shorter distance to target sequence
                //than our pattern. Now there r two options 1.) while converting pattern to target we would 'go through' the second barcode. In this
                //case minConversion/2 is always bigger than the editDist and we don t break
                //2.) the minconversion is the maximum conversion from pattern to second barcode and target is inbetween converting these two
                //now if the current barcode is however closer to target (minConversion/2), then this is the best match we cna ever find...
                if(bestEditDist < (minConversion/2))
                {
                    break;
                }
            }

        }

        //if we have several best matches we have to return false
        if(severalMatches)
        {
            return false;
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

        std::unordered_map<std::string, int> pattern_conversionrates;

        EdlibAlignConfig config;
};

class WildcardBarcode : public Barcode
{
    // TODO: move public parameter mismatches to a private and derived parameter
    //wildcardBarcode doe snot make use of mismatches yet, since anyways we do not know the sequence,
    //therefore its an unused parameter, just set for completeness as these classes derive from Barcode (initialized with mismatches, see up...)
    public:
    WildcardBarcode(int inMismatches, std::string name, int inLength) : Barcode(name, inMismatches, inLength) {}

    bool align(std::string& matchedBarcode, const std::string& target, const unsigned int positionInFastqLine,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)
    {
        (void)delNum; // silence unused parameter warning
        (void)insNum; // silence unused parameter warning
        (void)substNum; // silence unused parameter warning
        (void)reverse; // silence unused parameter warning

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

    bool align(std::string& matchedBarcode, const std::string& target, const unsigned int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)

        {
            (void)matchedBarcode; // silence unused parameter warning
            (void)target; // silence unused parameter warning
            (void)targetOffset; // silence unused parameter warning
            (void)targetEnd; // silence unused parameter warning
            (void)delNum; // silence unused parameter warning
            (void)insNum; // silence unused parameter warning
            (void)substNum; // silence unused parameter warning
            (void)reverse; // silence unused parameter warning

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

    bool align(std::string& matchedBarcode, const std::string& target, const unsigned int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)

        {
            (void)matchedBarcode; // silence unused parameter warning
            (void)target; // silence unused parameter warning
            (void)targetOffset; // silence unused parameter warning
            (void)targetEnd; // silence unused parameter warning
            (void)delNum; // silence unused parameter warning
            (void)insNum; // silence unused parameter warning
            (void)substNum; // silence unused parameter warning
            (void)reverse; // silence unused parameter warning

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

    bool align(std::string& matchedBarcode, const std::string& target, const unsigned int targetOffset,
        int& targetEnd, int& delNum, int& insNum, int& substNum,
        bool reverse = false)

        {
            (void)matchedBarcode; // silence unused parameter warning
            (void)target; // silence unused parameter warning
            (void)targetOffset; // silence unused parameter warning
            (void)targetEnd; // silence unused parameter warning
            (void)delNum; // silence unused parameter warning
            (void)insNum; // silence unused parameter warning
            (void)substNum; // silence unused parameter warning
            (void)reverse; // silence unused parameter warning

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