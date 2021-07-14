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
    virtual bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, fastqStats& stats) = 0;

};

class ConstantBarcode : public Barcode
{

    public:
    ConstantBarcode(std::string inPattern, int inMismatches) : pattern(inPattern),Barcode(inMismatches) {}
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, fastqStats& stats)
    {
        bool offset_shift = false;
        int tmpOffset = offset;
        if( !((offset-mismatches) < 0) )
        {
            tmpOffset = offset-mismatches;
            offset_shift = true;
        }
        int mismatch_length = (offset_shift ? 2*mismatches : mismatches);

        sequence = sequence.substr(tmpOffset, pattern.length() + mismatch_length);
        // e.g.: [AGTAGT]cccc: start=0 end=6 end is first not included idx
        if(levenshtein(sequence, pattern, mismatches, seq_start, seq_end, score))
        {
            if(offset_shift)
            {
                seq_start = seq_start-mismatches;
                seq_end = seq_end-mismatches;
            }
            return true;
        }
        else
        {
            //bigger mismatches
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
    bool match_pattern(std::string sequence, const int& offset, int& seq_start, int& seq_end, int& score, fastqStats& stats)
    {
        int match_count = 0;
        int best_start = 0, best_end = 0, best_score = mismatches+1;

        //calculate a new offset to include possible overlaps of barcodes in case of deletion and same
        //nucleotides at both ends
        bool offset_shift = false;
        int tmpOffset = offset;
        if( !((offset-mismatches) < 0) )
        {
            tmpOffset = offset-mismatches;
            offset_shift = true;
        }
        int mismatch_length = (offset_shift ? 2*mismatches : mismatches);

        for(std::string pattern : patterns)
        {
            int slice_end =  pattern.length() + mismatch_length;;
            std::string tmpSequence = sequence.substr(tmpOffset, slice_end);
            score = 0;
            seq_start = 0;
            seq_end = 0;
            if(levenshtein(tmpSequence, pattern, mismatches, seq_start, seq_end, score))
            {
                if(score < best_score)
                {
                    best_start = seq_start;
                    best_end = seq_end;
                    best_score = score;
                    match_count = 1;
                }
                else if(score == best_score)
                {
                    ++match_count;
                }
            }
        }

        if(match_count == 0)
        {
            return false;
        }
        else if(match_count > 1)
        {
            seq_start = best_start;
            seq_end = best_end;

            if(offset_shift)
            {
                seq_start = seq_start-mismatches;
                seq_end = seq_end-mismatches;
            }

            ++stats.multiBarcodeMatch;

            return false;
        }
        else
        {
            seq_start = best_start;
            seq_end = best_end;

            if(offset_shift)
            {
                seq_start = seq_start-mismatches;
                seq_end = seq_end-mismatches;
            }

            return true;
        }

    }

    private:
    std::vector<std::string> patterns;

};