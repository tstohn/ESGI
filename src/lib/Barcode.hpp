#pragma once

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
    //
    std::string generate_reverse_complement(std::string seq)
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
                default:
                    throw std::domain_error("Invalid nucleotide.");
            }
        };
        std::string newSeq = seq;
        std::transform(seq.crbegin(), seq.crend(), newSeq.begin(), lambda);
        return newSeq;
    }
    //overwritten function to match sequence pattern(s)
    virtual bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                            const int& differenceInBarcodeLength, bool startCorrection = false, bool reverse = false, bool fullLengthMapping = false) = 0;
    virtual std::vector<std::string> get_patterns() = 0;
    virtual bool is_wildcard() = 0;
    virtual bool is_constant() = 0;
};

class ConstantBarcode : public Barcode
{

    public:
    ConstantBarcode(std::string inPattern, int inMismatches) : pattern(inPattern),Barcode(inMismatches) 
    {
        revCompPattern = generate_reverse_complement(pattern);
    }
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       const int& differenceInBarcodeLength, bool startCorrection = false, bool reverse = false,
                       bool fullLengthMapping = false)
    {

        int tries = differenceInBarcodeLength;
        int tmpOffset = offset;
        bool offsetShiftBool = false;
        int offsetShiftValue = 0;

        //increment offset each round
        while(tries >= 0)
        {
            bool matchResult = private_match_pattern(sequence, tmpOffset, offsetShiftValue, seq_start, seq_end,score ,
                                     realBarcode, offsetShiftBool, startCorrection, reverse, fullLengthMapping);
            if(matchResult){return true;}

            //increment offset
            if( (++tmpOffset) > sequence.length() )
            {
                return false;
            }   
            offsetShiftBool = true;
            ++offsetShiftValue;
            --tries;
        }
        //if no mapping for each try worked
        return false;
    }
    std::vector<std::string> get_patterns()
    {
        std::vector<std::string> patterns = {pattern};
        return patterns;
    }
    bool is_wildcard(){return false;}
    bool is_constant(){return true;}

    private:
    bool private_match_pattern(std::string sequence, const int& offset, const int& offsetShiftValue, int& seq_start, int& seq_end, 
                               int& score, std::string& realBarcode, const bool& offsetShiftBool, bool startCorrection = false,
                               bool reverse = false, bool fullLengthMapping = false)
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
            //start and end are filled within levenshtein: if subsequences has a new offset this needs to be added to them
            if(offsetShiftBool)
            {
                seq_start = seq_start+offsetShiftValue;
                seq_end = seq_end+offsetShiftValue;
            }
            realBarcode = pattern;
            int diffEnd = pattern.length()-endInPattern;

            //calculate possible mapping that we oversaw bcs our window has a fixed size
            //if we have not mathced the whole pattern to subsequence and we can even still elongate to the end of the subsequence in the read
            if( (diffEnd != 0) && (sequence.length() >= offset + pattern.length() + diffEnd) )
            {
                std::string subSequenceElongated = sequence.substr(offset, pattern.length() + diffEnd);
                if(subSequence.length() != subSequenceElongated.length() && seq_end < subSequenceElongated.length())
                {
                    int extension = backBarcodeMappingExtension(subSequenceElongated, usedPattern, seq_end, endInPattern);
                    seq_end += extension;
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
};
class VariableBarcode : public Barcode
{

    public:
    VariableBarcode(std::vector<std::string> inPatterns, int inMismatches) : patterns(inPatterns), Barcode(inMismatches) 
    {
        for(std::string pattern : patterns)
        {
            std::string revCompPattern = generate_reverse_complement(pattern);
            revCompPatterns.push_back(revCompPattern);
        }
    }
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       const int& differenceInBarcodeLength, bool startCorrection = false,  bool reverse = false, bool fullLengthMapping = false)
    {
        int tries = differenceInBarcodeLength;
        int tmpOffset = offset;
        bool offsetShiftBool = false;
        int offsetShiftValue = 0;
        int tmpScore = mismatches + 1;
        int numberOfSameScoreResults = 0;
        while(tries >= 0)
        {
            int tmpNumberOfSameScoreResults = 0;
            bool matchResult = private_match_pattern(sequence, tmpOffset, offsetShiftValue, seq_start, seq_end,score ,
                                     realBarcode, offsetShiftBool, tmpNumberOfSameScoreResults, reverse, startCorrection);

            //if we already map sth to 100% don t bother and return true
            if(matchResult && score == 0){return true;}
            
            if( (++tmpOffset) > sequence.length() )
            {
                break;
            }
            offsetShiftBool = true;
            ++offsetShiftValue;
            --tries;

            //if we iterate through several windows, make sure we do not find multiple euqally good solutions
            //equally good solutions within ONE window are handled in private_match_pattern
            if(tmpScore == score)
            {
                ++numberOfSameScoreResults;
            }
            else if(tmpScore > score)
            {
                tmpScore = score;
                numberOfSameScoreResults = tmpNumberOfSameScoreResults;
            }

            //if no try results in a better mapping than mismatches plus 1, report false
            if( (tries < 0) && (tmpScore > mismatches))
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
        else if( (numberOfSameScoreResults == 0) && (tmpScore <= mismatches)){
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

    private:
    // IMPROVE FUNCTION:
    //      window size is: bases unmatched at end of previous sequence + mismatches for this barcode
    //      instead of one long seq with more allowed mismatches: move along a window and get barcode

    // first only for number of skipped bases call a sequences window; if this leads to nothing add more windows until number mismatches in barcode is reached
    bool private_match_pattern(std::string sequence, const int& offset, const int& offsetShiftValue, int& seq_start, int& seq_end, int& score, 
                               std::string& realBarcode, const bool& offsetShiftBool, int& numberOfSameScoreResults, 
                               bool reverse = false, bool startCorrection = false)
    {
        int match_count = 0;
        bool found_match = false;
        int best_start = 0, best_end = 0, best_score = mismatches+1;

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
                if(offsetShiftBool)
                {
                    seq_start = seq_start+offsetShiftValue;
                    seq_end = seq_end+offsetShiftValue;
                }
                int diffEnd = pattern.length()-endInPattern;
                //if we have not mathced the whole pattern to subsequence and we can even still elongate to the end of the subsequence in the read
                if( (diffEnd != 0) && (sequence.length() >= offset + pattern.length() + diffEnd) )
                {
                    std::string subSequenceElongated = sequence.substr(offset, pattern.length() + diffEnd);
                    if(subSequence.length() != subSequenceElongated.length() && seq_end < subSequenceElongated.length())
                    {
                        int extension = backBarcodeMappingExtension(subSequenceElongated, usedPattern, seq_end, endInPattern);
                        seq_end += extension;
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
            return false;
        }
        else if(match_count > 1)
        {
            seq_start = best_start;
            seq_end = best_end;
            score = best_score;

            ++numberOfSameScoreResults;
            return false;
        }
        else
        {
            seq_start = best_start;
            seq_end = best_end;
            score = best_score;

            numberOfSameScoreResults = 0;
            return true;
        }
    }
    std::vector<std::string> patterns;
    std::vector<std::string> revCompPatterns;

};

class WildcardBarcode : public Barcode
{
    // TODO: move public parameter mismatches to a private and derived parameter
    //wildcardBarcode doe snot make use of mismatches yet, since anyways we do not know the sequence,
    //therefore its an unused parameter, just set for completeness as these classes derive from Barcode (initialized with mismatches, see up...)
    public:
    WildcardBarcode(std::string inPattern, int inMismatches) : pattern(inPattern),Barcode(inMismatches) {}
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       const int& differenceInBarcodeLength, bool startCorrection = false, bool reverse = false, bool fullLengthMapping = false)
    {

        sequence = sequence.substr(offset, pattern.length());
        int end = (sequence.length() < pattern.length()) ? sequence.length() : pattern.length();
        // e.g.: [AGTAGT]cccc: start=0 end=6 end is first not included idx
        seq_start = 0;
        seq_end = end;
        realBarcode = sequence;
        return true;
    }
    std::vector<std::string> get_patterns()
    {
        std::vector<std::string> patterns = {pattern};
        return patterns;
    }
    bool is_wildcard(){return true;}
    bool is_constant(){return false;}

    private:
    std::string pattern; //just a string of "XXXXX"
};