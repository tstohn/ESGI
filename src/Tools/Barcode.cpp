#include "Barcode.hpp"

//Barcode::Barcode(int inMismatches) : mismatches(inMismatches) {}
//ConstantBarcode::ConstantBarcode(std::string inPattern, int inMismatches) : pattern(inPattern),Barcode(inMismatches) {}
//VariableBarcode::VariableBarcode(std::vector<std::string> inPatterns, int inMismatches) : patterns(inPatterns), Barcode(inMismatches) {}

/*
bool ConstantBarcode::match_pattern(std::string sequence, int offset, int seq_start, int seq_end, fastqStats& stats)
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

bool VariableBarcode::match_pattern(std::string sequence, int offset, int seq_start, int seq_end, fastqStats& stats)
{

    return true;
}*/