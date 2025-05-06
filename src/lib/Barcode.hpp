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
    virtual bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                               int& differenceInBarcodeLength, bool startCorrection = false, bool reverse = false, bool fullLengthMapping = false) = 0;
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
        std::string target = fastqLine.substr(targetOffset, pattern.length()+mismatches);

        //set the pattern to use for reverse or forward mapping
        std::string usedPattern = pattern;
        if(reverse){usedPattern = revCompPattern;}

        //map the pattern to the target sequence
        foundAlignment = run_alignment(usedPattern, target, targetEnd, config, delNum,  insNum, substNum);
        matchedBarcode = pattern;

        return foundAlignment;
    }

    /*
    ALGORITHM:
        iterate over patterns and match it to substring plus/minus mismatches on both sides
    allow mismatches at beginning and end (plus/minus those mismatches, bcs imagine a barcode match with a mismatch in the end,
    we then donnt know if its  really a msimatch, or a deletion of the barcode and already part of the next barcode...)
    each matched barcodes is described by the first and last match of the sequence (therefore can be shorter, than real sequence, but not longer)
    UMI or WildcardBarcodes are matched according to the two last matches in the neighboring sequences
    aligning by semi global alignment, bcs it could be that our pattern matches beyond the sequence, that is checked afterwards
    */
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       int& differenceInBarcodeLength, bool startCorrection = false, bool reverse = false,
                       bool fullLengthMapping = false)
    {

        int tries = differenceInBarcodeLength;
        int tmpOffset = offset;
        bool offsetShiftBool = false;
        int offsetShiftValue = 0;

        //increment offset each round
        bool matchFound = false;
        mappingSolution tmpSolution;

        //iterate over tries, and keep first found solution
        while(tries >= 0)
        {
            bool matchResult = private_match_pattern(sequence, tmpOffset, offsetShiftValue, seq_start, seq_end,score ,
                                     realBarcode, offsetShiftBool, differenceInBarcodeLength, startCorrection, reverse, fullLengthMapping);
            //if(matchResult){matchFound = true;return true;}
            
            //could be extended to iterate over all tries and only keep best solution
            if(matchResult && score < tmpSolution.score)
            {
                matchFound = true;
                tmpSolution = {seq_start, seq_end,score , realBarcode, differenceInBarcodeLength};
                break;
            }

            //increment offset
            if( (++tmpOffset) > sequence.length() && !matchFound )
            {
                return false;
            }   
            offsetShiftBool = true;
            ++offsetShiftValue;
            --tries;
        }

        //if no match found
        if(!matchFound){return false;}

        //else store best solution in origional variables
        seq_start = tmpSolution.seq_start;
        seq_end = tmpSolution.seq_end;
        score = tmpSolution.score;
        realBarcode = tmpSolution.realBarcode;
        differenceInBarcodeLength = tmpSolution.differenceInBarcodeLength;

        return true;
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
    bool private_match_pattern(std::string sequence, const int& offset, const int& offsetShiftValue, int& seq_start, int& seq_end, 
                               int& score, std::string& realBarcode, const bool& offsetShiftBool, int& diffEnd,
                               bool startCorrection = false, bool reverse = false, bool fullLengthMapping = false)
    {
        //set the pattern to use for reverse or forward mapping
        std::string usedPattern = pattern;
        if(reverse){usedPattern = revCompPattern;}

        std::string subSequence;
        if(!fullLengthMapping)
        {
            subSequence = sequence.substr(offset, pattern.length());
        }
        else
        {
            subSequence = sequence;
        }

        int endInPattern = 0; // store the number of missing bases in the pattern (in this case we might have to elongate the mapped sequence)
        // e.g.: [AGTAGT]cccc: start=0 end=6 end is first not included idx
        int startInPattern = 0;
        if(levenshtein(subSequence, usedPattern, mismatches, seq_start, seq_end, score, endInPattern, startInPattern))
        {
            //seq_start starts potentially with 1, and seq_end is in a perfect match length (zero row, col is filled with zeroes in edit-dist)
            int differencePatternLengthMappingLength = pattern.length()-(seq_end);
            diffEnd = pattern.length()-endInPattern;
            //we only set the offset based on mapping differencen when: 
            //1.) the sequence is not mapped fully to the pattern (therefore, there must have been deletion) 
            //2.) the pattern is mapped till the end (below)
            if(diffEnd!= 0){differencePatternLengthMappingLength = 0;}

            //start and end are filled within levenshtein: if subsequences has a new offset this needs to be added to them
            if(offsetShiftBool)
            {
                seq_start = seq_start+offsetShiftValue;
                seq_end = seq_end+offsetShiftValue;
            }

            realBarcode = pattern;

            //calculate possible mapping that we oversaw bcs our window has a fixed size
            //if we have not mathced the whole pattern to subsequence and we can even still elongate to the end of the subsequence in the read
            if( (differencePatternLengthMappingLength == 0) && (diffEnd > 0) && (sequence.length() >= offset + pattern.length() + differencePatternLengthMappingLength) )
            {
                std::string subSequenceElongated = sequence.substr(offset, pattern.length() + diffEnd);
                if(subSequence.length() != subSequenceElongated.length() && seq_end < subSequenceElongated.length())
                {
                    int extension = backBarcodeMappingExtension(subSequenceElongated, usedPattern, seq_end, endInPattern);
                    seq_end += extension;
                    diffEnd -= extension;
                }
            }

            //possible extension at the beginning: is only possible if before we had a wildcard
            if(startInPattern > 0 && startCorrection && (offset >= startInPattern) )
            {
                //we have startInPattern additional bases to check before sequence
                std::string subSequenceElongated = sequence.substr(offset-startInPattern, pattern.length());
                int extension = frontBarcodeMappingExtension(subSequenceElongated, usedPattern, seq_start, startInPattern);
                seq_start -= extension;
            }

            return true;
        }
        else
        {
            //bigger mismatches
            return false;
        }
    }
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
        if(equalLengthBarcodes)
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

            delNum=insNum=substNum=targetEnd=0;
            //get pattern, its length can vary
            std::string usedPattern = patternsToMap.at(patternIdx);

            //define target sequence (can differ for every barcode due to its length)
            std::string target = fastqLine.substr(targetOffset, usedPattern.length()+mismatches);

            //std::cout << " in barcode trying barcode: " << usedPattern << " with target sequ " << target<< "\n";

            //map the pattern to the target sequence
            foundAlignment = run_alignment(usedPattern, target, targetEnd, config, delNum,  insNum, substNum);
            
            if(foundAlignment && (delNum+insNum+substNum)<bestEditDist)
            {
                bestFoundAlignment = foundAlignment;
                bestFoundPattern = patterns.at(patternIdx); //the barcode is the TRUE forward barcode, not the reverse complement
                bestTargetEnd = targetEnd;
                bestEditDist = (delNum+insNum+substNum);
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

    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       int& differenceInBarcodeLength, bool startCorrection = false,  bool reverse = false, bool fullLengthMapping = false)
    {

        int tries = differenceInBarcodeLength;
        int tmpOffset = offset;
        bool offsetShiftBool = false;
        int offsetShiftValue = 0;
        int numberOfSameScoreResults = 0;
        mappingSolution tmpSolution;
        //iterate over tries, and keep first found solution
        while(tries >= 0)
        {
            int tmpNumberOfSameScoreResults = 0;
            bool matchResult = private_match_pattern(sequence, tmpOffset, offsetShiftValue, seq_start, seq_end,score ,
                                     realBarcode, offsetShiftBool, tmpNumberOfSameScoreResults, differenceInBarcodeLength, reverse, startCorrection);

            //could be extended to iterate over all tries and only keep best solution
            if(matchResult && score < tmpSolution.score)
            {
                tmpSolution = {seq_start, seq_end,score , realBarcode, differenceInBarcodeLength};
                numberOfSameScoreResults = tmpNumberOfSameScoreResults;

                break;
            }
            //if we iterate through several windows, make sure we do not find multiple euqally good solutions
            //equally good solutions within ONE window are handled in private_match_pattern
            else if( matchResult && (tmpSolution.score == score) && (tmpSolution.realBarcode != realBarcode))
            {
                ++numberOfSameScoreResults;
            }
            
            if( (++tmpOffset) > sequence.length() )
            {
                break;
            }
            offsetShiftBool = true;
            ++offsetShiftValue;
            --tries;

            //if no try results in a better mapping than mismatches plus 1, report false
            if( (tries < 0) && (tmpSolution.score > mismatches))
            {
                return false;
            } //if also the last try failed
        }
        //if the best match has several times this score, report also false and a MultiBarcodeMatch
        if(numberOfSameScoreResults > 0)
        {            
            //we mapped several barcodes with the same score
            return false;
        }
        else if( (numberOfSameScoreResults == 0) && (tmpSolution.score <= mismatches))
        {

            //else store best solution in origional variables
            seq_start = tmpSolution.seq_start;
            seq_end = tmpSolution.seq_end;
            score = tmpSolution.score;
            realBarcode = tmpSolution.realBarcode;
            differenceInBarcodeLength = tmpSolution.differenceInBarcodeLength;

            return true;
        }

        return false;
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
    // IMPROVE FUNCTION:
    //      window size is: bases unmatched at end of previous sequence + mismatches for this barcode
    //      instead of one long seq with more allowed mismatches: move along a window and get barcode

    // first only for number of skipped bases call a sequences window; if this leads to nothing add more windows until number mismatches in barcode is reached
    bool private_match_pattern(std::string sequence, const int& offset, const int& offsetShiftValue, int& seq_start, int& seq_end, int& score, 
                               std::string& realBarcode, const bool& offsetShiftBool, int& numberOfSameScoreResults, int& diffEnd,
                               bool reverse = false, bool startCorrection = false)
    {

        int match_count = 0;
        bool found_match = false;
        int best_start = 0, best_end = 0, best_score = mismatches+1, tmpDiff = 0, bestDiff = 0;

        std::vector<std::string> patternsToMap = patterns;
        if(reverse){patternsToMap = revCompPatterns;}

        for(int patternIdx = 0; patternIdx!= patternsToMap.size(); ++patternIdx)
        {
            std::string pattern = patterns.at(patternIdx);
            std::string usedPattern = patternsToMap.at(patternIdx);

            std::string subSequence = sequence.substr(offset, pattern.length());

            score = 0;
            seq_start = 0;
            seq_end = 0;
            int endInPattern = 0; // store the number of missing bases in the pattern (in this case we might have to elongate the mapped sequence)
            int startInPattern = 0;

            if(levenshtein(subSequence, usedPattern, mismatches, seq_start, seq_end, score, endInPattern, startInPattern))
            {
                //seq_start starts potentially with 1, and seq_end is in a perfect match length (edit dist has zeroes in first row, col)
                int differencePatternLengthMappingLength = pattern.length()-(seq_end);
                tmpDiff = pattern.length()-endInPattern;
                //we only set the offset based on mapping differencen when: 
                //1.) the sequence is not mapped fully to the pattern (therefore, there must have been deletion) 
                //2.) the pattern is mapped till the end (below)
                if(tmpDiff!= 0){differencePatternLengthMappingLength = 0;}

                if(offsetShiftBool)
                {
                    seq_start = seq_start+offsetShiftValue;
                    seq_end = seq_end+offsetShiftValue;
                }

                //int differencePatternLengthMappingLength = pattern.length() - seq_end;
                //if we have not mathced the whole pattern to subsequence and we can even still elongate to the end of the subsequence in the read
                if( (differencePatternLengthMappingLength == 0) && (tmpDiff > 0) && (sequence.length() >= offset + pattern.length() + differencePatternLengthMappingLength) )
                {
                    std::string subSequenceElongated = sequence.substr(offset, pattern.length() + tmpDiff);
                    if(subSequence.length() != subSequenceElongated.length() && seq_end < subSequenceElongated.length())
                    {
                        int extension = backBarcodeMappingExtension(subSequenceElongated, usedPattern, seq_end, endInPattern);
                        seq_end += extension;
                        tmpDiff -= extension;
                    }
                }

                if( (startInPattern > 0) && (offset >=startInPattern) && startCorrection )
                {
                    //we have startInPattern additional bases to check before sequence: is only possible if before we had a wildcard
                    std::string subSequenceElongated = sequence.substr(offset-startInPattern, pattern.length());

                    int extension = frontBarcodeMappingExtension(subSequenceElongated, usedPattern, seq_start, startInPattern);
                    seq_start -= extension;
                }
                //score comparison
                if( (score < best_score))
                {
                    best_start = seq_start;
                    best_end = seq_end;
                    best_score = score;
                    realBarcode = pattern;
                    match_count = 1;
                    found_match = true;
                    bestDiff = tmpDiff;
                }
                else if( (score == best_score) & found_match)
                {
                    ++match_count;
                }
            }

        }
        //compare all results for different patterns
        if(match_count == 0)
        {
            numberOfSameScoreResults = 0;
            diffEnd = 0;
            return false;
        }
        else if(match_count > 1)
        {
            seq_start = best_start;
            seq_end = best_end;
            score = best_score;
            diffEnd = bestDiff;

            ++numberOfSameScoreResults;
            return false;
        }
        else
        {
            seq_start = best_start;
            seq_end = best_end;
            score = best_score;
            diffEnd = bestDiff;

            numberOfSameScoreResults = 0;
            return true;
        }
    }
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
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       int& differenceInBarcodeLength, bool startCorrection = false, bool reverse = false, bool fullLengthMapping = false)
    {

        sequence = sequence.substr(offset, length);
        int end = (sequence.length() < length) ? sequence.length() : length;
        // e.g.: [AGTAGT]cccc: start=0 end=6 end is first not included idx
        seq_start = 0;
        seq_end = end;
        realBarcode = sequence;
        return true;
    }
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
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       int& differenceInBarcodeLength, bool startCorrection = false, bool reverse = false, bool fullLengthMapping = false)
    {
        return false;
    }
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
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       int& differenceInBarcodeLength, bool startCorrection = false, bool reverse = false, bool fullLengthMapping = false)
    {
        return false;
    }
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
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       int& differenceInBarcodeLength, bool startCorrection = false, bool reverse = false, bool fullLengthMapping = false)
    {
        return false;
    }
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