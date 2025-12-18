#pragma once

#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

#include "Barcode.hpp"

//a class for temporary results of a single fastq-line
//this information is then merged into the final statistic object
struct OneLineDemultiplexingStats
{

        //for every pattern it saves until where we could map, only in case we could not map a line
        //do we transfer this information into the final structure;
        //string == patternName, int = position within this pattern which was not mapped and where mapping failed
        //it only saves the positions for failed lines or lines where we did not map all patterns (no missing patterns
        //between fw and rv reads, we might still map perfectly if the overlap is a constant barcode)
        //if the last few bases in a read are also shorter than the next barcode that has to be mapped we assign this read to be
        //perfectly mapped (the seqToShort check does therefore not lead to failed reads)
        //a simple integer would also be enough at the moment since we only allow this quality feature when a single pattern is provided
        //but for potential extension later on we keep it a pair that aslso store the actuakl barcode it comes from
        std::pair<std::string, int> failedLinesMappingFw;
        std::pair<std::string, int> failedLinesMappingRv;

        //GLOBAL PARAMETERS
        //parameters that are evaluated over the whole fastq line
        //e.g. perfect match occurs only if ALL barcodes match perfectly in a fastq line
        unsigned long long perfectMatches = 0;
        unsigned long long moderateMatches = 0;

        //MISMATCH TYPE PARAMETERS
        //simply mismatches in order of positions, we fill the final dict after every read alignment
        //IMPORTANT: fill this vector for every barcode positions (also UMIs and leave it at zero), 
        //we later 'pick' only the valid barcode positions
        //e.g.: [3,0,0,5,1] for a pattern like [BC1,UMI,DNA,LINKER,BC2]
        //the last mismatches are always substitutions in case it can also be deletions
        //in case the FW and RV reads do not overlap, and the missing sequence is only constant barcodes, we 
        //do not count the MM in this barcode
        std::vector<int> insertions;
        std::vector<int> deletions;
        std::vector<int> substitutions;

};
typedef std::shared_ptr<OneLineDemultiplexingStats> OneLineDemultiplexingStatsPtr;

class DemultiplexingStats{

    public:     

        void initializeStats(const MultipleBarcodePatternVectorPtr& barcodePatternList);
        void update(OneLineDemultiplexingStatsPtr lineStatsPtr, bool result, std::string& foundPatternName, std::vector<std::string>& barcodeList);
        
        void write_mm_number(const std::string& outputFile);
        void write_failing_position(const std::string& outputFile);
        void write_mm_types(const std::string& outputFile);
        void write(const std::string& directory, const std::string& prefix, const int patternNumber);

        //UPDATE FUNCTIONS
        void update_failedLinesMapping(const std::pair<std::string, int>& failedFw, const std::pair<std::string, int>& failedRv)
        {

            //if there was an error
            //update the last mapped position on FW read
            if(failedFw.first != "")
            {
                //key is <PATTERN>_<LastMappedPosition>
                std::string keyFw = failedFw.first + "_" + std::to_string(failedFw.second);
                ++failedLinesMappingFw.at(keyFw);
            }

            //if there was an error
            //update the last mapped position on RV read
            if(failedRv.first != "")
            {
                //key is <PATTERN>_<LastMappedPosition>
                std::string keyRv = failedRv.first + "_" + std::to_string(failedRv.second);
                ++failedLinesMappingRv.at(keyRv);
            }
        }

        //update mapping statistics: perfect/ moderate mapping is only set when mapped
        //if not mapped both are zero
        void update_global_parameters(bool result, OneLineDemultiplexingStatsPtr lineStatsPtr)
        {
            if(!result)
            {
                ++noMatches;
            }
            else
            {
                perfectMatches += lineStatsPtr->perfectMatches;
                moderateMatches += lineStatsPtr->moderateMatches;
            }
        }

        void update_mismatch_types(OneLineDemultiplexingStatsPtr lineStatsPtr, std::string& foundPatternName, std::vector<std::string>& barcodeList)
        {
            //key is a combinations of <PATTERN>_<POSITION>_<BARCODE>
            // temporary input (lineStatsPtr->insertions) is a vector of e.g., insertions at every barcode position
            //update insertions
            for (int validBarcodePos : validPositions.at(foundPatternName)) 
            {
                std::string key = foundPatternName + "_" + std::to_string(validBarcodePos) + "_" + barcodeList.at(validBarcodePos);
                insertions.at(key) += lineStatsPtr->insertions.at(validBarcodePos);
            }

            //update deletions
            for (int validBarcodePos : validPositions.at(foundPatternName)) 
            {
                std::string key = foundPatternName + "_" + std::to_string(validBarcodePos) + "_" + barcodeList.at(validBarcodePos);
                deletions.at(key) += lineStatsPtr->deletions.at(validBarcodePos);
            }

            //update substitutions
            for (int validBarcodePos : validPositions.at(foundPatternName)) 
            {
                std::string key = foundPatternName + "_" + std::to_string(validBarcodePos) + "_" + barcodeList.at(validBarcodePos);
                substitutions.at(key) += lineStatsPtr->substitutions.at(validBarcodePos);
            }

        }

        void update_mismatch_numbers(OneLineDemultiplexingStatsPtr lineStatsPtr, std::string& foundPatternName, std::vector<std::string>& barcodeList)
        {
            //key is a combinations of <PATTERN>_<POSITION>_<BARCODE>

            //update how often we see which number of MM for a key
            //only use valid barcode positions (no UMI, DNA, STOP etc positions)
            for (int validBarcodePos : validPositions.at(foundPatternName)) 
            {
                //generate key
                std::string key = foundPatternName + "_" + std::to_string(validBarcodePos) + "_" + barcodeList.at(validBarcodePos);
                                
                //sum up insertions/ deletions/ substitutions at every positions
                int mmSum = lineStatsPtr->insertions.at(validBarcodePos) + 
                            lineStatsPtr->deletions.at(validBarcodePos) + 
                            lineStatsPtr->substitutions.at(validBarcodePos);

                //the value of the map mismatchNumber is a vector of the length of potential mismatches +1
                //we need to increment the count at the specific position for the number of mismatches
                ++(mismatchNumber.at(key)).at(mmSum);

            }
        }

        void combine_statistics(std::vector<std::shared_ptr<DemultiplexingStats>>& statisticsList);

        //GETTER FUNCTIONS
        ///number of perfect matches
        unsigned long long get_perfect_matches() const
        {
            return perfectMatches;
        }
        ///number of matches with mismatches
        unsigned long long get_moderat_matches() const
        {
            return moderateMatches;
        }
        ///number of failed matches
        unsigned long long get_failed_matches() const
        {
            return noMatches;
        }
        ///the dictionary of mismatches per barcode
        const std::map<std::string, std::vector<int> > get_mismatch_dict()
        {
            return mismatchNumber;
        }

    private:

        //adding quality scores to class
        //this is only stored for failed lines:
        //imagine a line can not be mapped to pattern A
        //then it stores how far we could map in pattern A 
        //key is <PATTERN>_<POSITION>
        //so it might be extended in the future to store this for all patterns, for only it only works when we have ONE pattern
        std::map<std::string, int> failedLinesMappingFw;
        std::map<std::string, int> failedLinesMappingRv;

        //GLOBAL PARAMETERS
        //parameters that are evaluated over the whole fastq line
        //e.g. perfect match occurs only if ALL barcodes match perfectly in a fastq line
        unsigned long long perfectMatches = 0;
        unsigned long long noMatches = 0;
        unsigned long long moderateMatches = 0;

        //not for all barcodes can we store mismatches, etc.
        //therefore we need to keep track of the barcode positions that can actually contain MM (exclude all UMI, STOP, DNA, ETC. pos)
        //key <PATTERN> -> value vector of valid positions
        //the valid barcode positions are positions within the vector of barcode (so already excluding all stop, dna, read end barcodes)
        std::map<std::string, std::vector<int>> validPositions;
        
        //MISMATCH TYPE PARAMETERS
        //key is a combinations of <PATTERN>_<POSITION>_<BARCODE>
        std::map<std::string, int> insertions;
        std::map<std::string, int> deletions;
        std::map<std::string, int> substitutions;

        //MISMATCH NUMBER PARAMETERS
        //key is a combinations of <PATTERN>_<POSITION>_<BARCODE>
        //<key> : [0MM, 1MM, 2MM, ..., MAX_MM]
        std::map<std::string, std::vector<int>> mismatchNumber;

};
