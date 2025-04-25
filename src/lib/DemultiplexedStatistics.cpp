#include "DemultiplexedStatistics.hpp"

//initialize the stats dictionary of mismatches per barcode
//initialize the dictionary of mismatches for each possible barcode
// the vector in the dict has length "mismatches + 2" for each barcode
// one entry for zero mismatches, eveery number from, 1 to mismatches and 
//one for more mismatches than the max allowed number of mismathces
void DemultiplexingStats::initializeStats(const MultipleBarcodePatternVectorPtr& barcodePatternList)
{
    //initialize lock
    statsLock = std::make_unique<std::mutex>();

    //INITIALiZE THE FAILED LINE DICTIONARIES (one for FW and RV): <PATTERN_NAME> => <LAST_POSITION_MAPPED>
    //for pattern
    for(BarcodePatternPtr patternPtr : *barcodePatternList)
    {
        //for barcode within this pattern
        for(int barcodePos = 0; barcodePos < patternPtr->barcodePattern->size(); ++barcodePos)
        {
            //initialize the failed line positions for fW and rv reads
            std::string pattern_position = patternPtr->patternName + "_" + std::to_string(barcodePos);
            failedLinesMappingFw.insert(std::make_pair(pattern_position, 0));
            failedLinesMappingRv.insert(std::make_pair(pattern_position, 0));

            int mismatches = patternPtr->barcodePattern->at(barcodePos)->mismatches;
            const std::vector<std::string> variableBarcodesVec = patternPtr->barcodePattern->at(barcodePos)->get_patterns();

            if(!patternPtr->barcodePattern->at(barcodePos)->is_wildcard() && 
               !patternPtr->barcodePattern->at(barcodePos)->is_stop() && 
               !patternPtr->barcodePattern->at(barcodePos)->is_dna() && 
               !patternPtr->barcodePattern->at(barcodePos)->is_read_end())
                {
                    //this is a valid position in this pattern, safe it so we can fill it later on
                    validPositions[patternPtr->patternName].push_back(barcodePos);

                    //for variable barcodes // (or single constant one)
                    for(const std::string& barcodeString : variableBarcodesVec)
                    {
                        //create key for dictionaries: <pattern=name>_<barcodePosition>_<actual Barcode>
                        std::string pattern_position_barcode = patternPtr->patternName + "_" + std::to_string(barcodePos) + "_" + barcodeString;

                        //INITIALiZE THE MISMATCH TYPE DICTIONARY
                        insertions.insert(std::make_pair(pattern_position_barcode, 0));
                        deletions.insert(std::make_pair(pattern_position_barcode, 0));
                        substitutions.insert(std::make_pair(pattern_position_barcode, 0));

                        //INITIALiZE THE MISMATCH NUMBER DICTIONARY
                        std::vector<int> mismatchVector(mismatches + 1, 0);
                        mismatchNumber.insert(std::make_pair(pattern_position_barcode, mismatchVector));
                    }
                }
        }   
    }
}  

void DemultiplexingStats::update(OneLineDemultiplexingStatsPtr lineStatsPtr, bool result, std::string& foundPatternName, std::vector<std::string>& barcodeList)
{

    //lock statistics and update, this object is updated across all threads in one shared instance
    std::lock_guard<std::mutex> guard(*statsLock);
    //update weather the line was mapped perfectly, moderately, or not at all
    update_global_parameters(result, lineStatsPtr, foundPatternName);
    
    //result dependent updates:
    // 1.) for mapped lines the mismatch types and numbers
    // 2.) for failed lines until where we mapped
    if(result)
    {
        //1.)
        update_mismatch_types(lineStatsPtr, foundPatternName, barcodeList);
        update_mismatch_numbers(lineStatsPtr, foundPatternName, barcodeList);
    }
    else
    {
        //2.)
        //update for every pattern until where we could map (last mapped barcode position)
        update_failedLinesMapping(lineStatsPtr->failedLinesMappingFw, lineStatsPtr->failedLinesMappingRv);
    }

}

void DemultiplexingStats::write(const std::string& directory)
{
    //remove files if they already exist
    std::string barcodeMismatchNumber = directory + "/Quality_numberMM.txt";
    std::string barcodeMismatchType = directory + "/Quality_typeMM.txt";
    std::string barcodeLastPosMapped = directory + "/Quality_lastPositionMapped.txt";

    // remove outputfile if it exists
    std::remove(barcodeMismatchNumber.c_str());
    std::remove(barcodeMismatchType.c_str());
    std::remove(barcodeLastPosMapped.c_str());

    //fill the number of MM

    //fill the type of MM

    //fill last mapped positions

}