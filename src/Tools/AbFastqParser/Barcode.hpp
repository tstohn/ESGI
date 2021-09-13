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
    //overwritten function to match sequence pattern(s)
    virtual bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                               const int& differenceInBarcodeLength, fastqStats& stats, bool startCorrection = false) = 0;
    virtual std::vector<std::string> get_patterns() = 0;
    virtual bool is_wildcard() = 0;

};

class ConstantBarcode : public Barcode
{

    public:
    ConstantBarcode(std::string inPattern, int inMismatches) : pattern(inPattern),Barcode(inMismatches) {}
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       const int& differenceInBarcodeLength, fastqStats& stats, bool startCorrection = false)
    {
        //if(!private_match_pattern(sequence, offset, seq_start, seq_end,score ,realBarcode ,stats, false))
        //{
         //   return(private_match_pattern(sequence, offset, seq_start, seq_end,score ,realBarcode ,stats, true));
        //}

        int tries = differenceInBarcodeLength;
        int tmpOffset = offset;
        bool offsetShiftBool = false;
        int offsetShiftValue = 0;
        //increment offset each round
        while(tries >= 0)
        {
            bool matchResult = private_match_pattern(sequence, tmpOffset, offsetShiftValue, seq_start, seq_end,score ,
                                     realBarcode ,stats, offsetShiftBool, startCorrection);
            if(matchResult){return true;}

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

    private:
    bool private_match_pattern(std::string sequence, const int& offset, const int& offsetShiftValue, int& seq_start, int& seq_end, 
                               int& score, std::string& realBarcode, fastqStats& stats, const bool& offsetShiftBool, bool startCorrection = false)
    {
     
        std::string subSequence = sequence.substr(offset, pattern.length());

        int endInPattern = 0; // store the number of missing bases in the pattern (in this case we might have to elongate the mapped sequence)
        // e.g.: [AGTAGT]cccc: start=0 end=6 end is first not included idx
        int startInPattern = 0;
        if(levenshtein(subSequence, pattern, mismatches, seq_start, seq_end, score, endInPattern, startInPattern))
        {
            if(offsetShiftBool)
            {
                seq_start = seq_start+offsetShiftValue;
                seq_end = seq_end+offsetShiftValue;
            }
            realBarcode = pattern;
            int diffEnd = pattern.length()-endInPattern;
            if(diffEnd != 0)
            {
                std::string subSequenceElongated = sequence.substr(offset, pattern.length() + diffEnd);
                if(subSequence.length() != subSequenceElongated.length() && seq_end < subSequenceElongated.length())
                {
                    int extension = backBarcodeMappingExtension(subSequenceElongated, pattern, seq_end, endInPattern);
                    seq_end += extension;
                }
            }
            if(startInPattern > 0 && startCorrection)
            {
                //we have startInPattern additional bases to check before sequence
                std::string subSequenceElongated = sequence.substr(offset-startInPattern, pattern.length());
                int extension = frontBarcodeMappingExtension(subSequenceElongated, pattern, seq_start, startInPattern);
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
};
class VariableBarcode : public Barcode
{

    public:
    VariableBarcode(std::vector<std::string> inPatterns, int inMismatches) : patterns(inPatterns), Barcode(inMismatches) {}
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       const int& differenceInBarcodeLength, fastqStats& stats, bool startCorrection = false)
    {
        int tries = differenceInBarcodeLength;
        int tmpOffset = offset;
        bool offsetShiftBool = false;
        int offsetShiftValue = 0;
        int tmpScore = mismatches + 1;
        int numberOfSameScoreResults = 0;
        while(tries >= 0)
        {
            bool matchResult = private_match_pattern(sequence, tmpOffset, offsetShiftValue, seq_start, seq_end,score ,
                                     realBarcode ,stats, offsetShiftBool, numberOfSameScoreResults, startCorrection);

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
                tmpScore = score;

            if( (tries < 0) && (score > mismatches)){return false;} //if also the last try failed
        }

        if( (numberOfSameScoreResults == 0) && (score <= mismatches)){return true;}
        //if a match pattern worked return true
        return false;
    }
    std::vector<std::string> get_patterns()
    {
        return patterns;
    }
    bool is_wildcard(){return false;}

    private:
    // IMPROVE FUNCTION:
    //      window size is: bases unmatched at end of previous sequence + mismatches for this barcode
    //      instead of one long seq with more allowed mismatches: move along a window and get barcode

    // first only for number of skipped bases call a sequences window; if this leads to nothing add more windows until number mismatches in barcode is reached
    bool private_match_pattern(std::string sequence, const int& offset, const int& offsetShiftValue, int& seq_start, int& seq_end, int& score, 
                               std::string& realBarcode, fastqStats& stats, const bool& offsetShiftBool, int& numberOfSameScoreResults, bool startCorrection = false)
    {
        int match_count = 0;
        bool found_match = false;
        int best_start = 0, best_end = 0, best_score = mismatches+1;

        for(std::string pattern : patterns)
        {
            std::string subSequence = sequence.substr(offset, pattern.length());

            score = 0;
            seq_start = 0;
            seq_end = 0;
            int endInPattern = 0; // store the number of missing bases in the pattern (in this case we might have to elongate the mapped sequence)
            //std::cout << "# " << pattern << "\n";
            int startInPattern = 0;

            if(levenshtein(subSequence, pattern, mismatches, seq_start, seq_end, score, endInPattern, startInPattern))
            {
                if(offsetShiftBool)
                {
                    seq_start = seq_start+offsetShiftValue;
                    seq_end = seq_end+offsetShiftValue;
                }
                int diffEnd = pattern.length()-endInPattern;
                if(diffEnd != 0)
                {
                    std::string subSequenceElongated = sequence.substr(offset, pattern.length() + diffEnd);
                    if(subSequence.length() != subSequenceElongated.length() && seq_end < subSequenceElongated.length())
                    {
                        int extension = backBarcodeMappingExtension(subSequenceElongated, pattern, seq_end, endInPattern);
                        seq_end += extension;
                    }
                }
                if( (startInPattern > 0) && (offset >startInPattern) && startCorrection )
                {
                    //we have startInPattern additional bases to check before sequence
                    std::string subSequenceElongated = sequence.substr(offset-startInPattern, pattern.length());

                    int extension = frontBarcodeMappingExtension(subSequenceElongated, pattern, seq_start, startInPattern);
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

            if(offsetShiftBool)
            {
                seq_start = seq_start-mismatches;
                seq_end = seq_end-mismatches;
            }

            ++stats.multiBarcodeMatch;
            ++numberOfSameScoreResults;

            return false;
        }
        else
        {
            seq_start = best_start;
            seq_end = best_end;
            score = best_score;

            if(offsetShiftBool)
            {
                seq_start = seq_start-mismatches;
                seq_end = seq_end-mismatches;
            }
            numberOfSameScoreResults = 0;

            return true;
        }
    }
    std::vector<std::string> patterns;

};

class WildcardBarcode : public Barcode
{
    // TODO: move public parameter mismatches to a private and derived parameter
    //wildcardBarcode doe snot make use of mismatches yet, since anyways we do not know the sequence,
    //therefore its an unused parameter, just set for completeness as these classes derive from Barcode (initialized with mismatches, see up...)
    public:
    WildcardBarcode(std::string inPattern, int inMismatches) : pattern(inPattern),Barcode(inMismatches) {}
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, std::string& realBarcode, 
                       const int& differenceInBarcodeLength, fastqStats& stats, bool startCorrection = false)
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

    private:
    std::string pattern;
};