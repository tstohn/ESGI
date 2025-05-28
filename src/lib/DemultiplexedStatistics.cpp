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
        int actualPatternPos = 0; //when demultiplexing the line we only push certain barcodes inot the demultiplexed list
        //we DO NOT push stop/ read-end or DNA patterns in there, and therefore also CAN NOT count these indices as indices
        //that store patterns with potential DEL/INS/MIS

        for(size_t barcodePos = 0; barcodePos < patternPtr->barcodePattern->size(); ++barcodePos)
        {

            //initialize the failed line positions for fW and rv reads
            std::string pattern_position = patternPtr->patternName + "_" + std::to_string(barcodePos);
            failedLinesMappingFw.insert(std::make_pair(pattern_position, 0));
            failedLinesMappingRv.insert(std::make_pair(pattern_position, 0));

            int mismatches = patternPtr->barcodePattern->at(barcodePos)->mismatches;
            const std::vector<std::string> variableBarcodesVec = patternPtr->barcodePattern->at(barcodePos)->get_patterns();

            //if pattern is a wildcard, then this pattern WILL BE in the demultiplexed barcode list, therefore the actualPatternPos
            //has to increase, however, we do not store this as a valid abrcode position with IN/DEL/MIS
            if(patternPtr->barcodePattern->at(barcodePos)->is_wildcard()){++actualPatternPos;}
            else if(!patternPtr->barcodePattern->at(barcodePos)->is_stop() && 
               !patternPtr->barcodePattern->at(barcodePos)->is_dna() && 
               !patternPtr->barcodePattern->at(barcodePos)->is_read_end())
                {

                    //this is a valid position in this pattern, safe it so we can fill it later on
                    validPositions[patternPtr->patternName].push_back(actualPatternPos);

                    //for variable barcodes // (or single constant one)
                    for(const std::string& barcodeString : variableBarcodesVec)
                    {
                        //create key for dictionaries: <pattern=name>_<barcodePosition>_<actual Barcode>
                        std::string pattern_position_barcode = patternPtr->patternName + "_" + std::to_string(actualPatternPos) + "_" + barcodeString;

                        //INITIALiZE THE MISMATCH TYPE DICTIONARY
                        insertions.insert(std::make_pair(pattern_position_barcode, 0));
                        deletions.insert(std::make_pair(pattern_position_barcode, 0));
                        substitutions.insert(std::make_pair(pattern_position_barcode, 0));

                        //INITIALiZE THE MISMATCH NUMBER DICTIONARY
                        std::vector<int> mismatchVector(mismatches + 1, 0);
                        mismatchNumber.insert(std::make_pair(pattern_position_barcode, mismatchVector));
                    }

                    ++actualPatternPos;
                }
        }   
    }
}  

void DemultiplexingStats::update(OneLineDemultiplexingStatsPtr lineStatsPtr, bool result, std::string& foundPatternName, std::vector<std::string>& barcodeList)
{

    //lock statistics and update, this object is updated across all threads in one shared instance
    std::lock_guard<std::mutex> guard(*statsLock);
    //update weather the line was mapped perfectly, moderately, or not at all
    update_global_parameters(result, lineStatsPtr);

    //result dependent updates:
    // 1.) for mapped lines the mismatch types and numbers
    // 2.) for failed lines until where we mapped
    if(result)
    {
        //1.)
        update_mismatch_types(lineStatsPtr, foundPatternName, barcodeList);
        update_mismatch_numbers(lineStatsPtr, foundPatternName, barcodeList);
    }
    else if(lineStatsPtr != nullptr) //it is a nullptr in case we have several patterns
    {
        //2.)
        //update for every pattern until where we could map (last mapped barcode position)
        //this can only be updated if we have only one pattern, otherwise we would have to create this per pattern
        update_failedLinesMapping(lineStatsPtr->failedLinesMappingFw, lineStatsPtr->failedLinesMappingRv);
    }

}

void DemultiplexingStats::write_mm_types(const std::string& outputFile) 
{
    std::ofstream out(outputFile);
    if (!out) {
        std::cerr << "Could not open file for writing number of insertions: " << outputFile << "\n";
        return;
    }

    // Write header
    out << "PATTERN\tPOSITION\tBARCODE\tMISMATCH_TYPE\tCOUNT\n";

    //write insertions
    for (const auto& [key, value] : insertions) {
        std::stringstream ss(key);
        std::string segment;
        std::vector<std::string> parts;

        bool elementIsInBrackets = false;
        std::string bracketString = "";
        while (std::getline(ss, segment, '_')) 
        {
            if(!elementIsInBrackets && segment.at(0) == '\'')
            {
                if(segment.back() == '\'')
                {
                    parts.push_back(segment);
                }
                else
                {
                    elementIsInBrackets = true;
                    bracketString += segment + "_";
                }
            }
            else if(elementIsInBrackets && segment.back() == '\'')
            {
                assert(elementIsInBrackets);
                elementIsInBrackets = false;
                bracketString += segment;
                parts.push_back(bracketString);
                bracketString = "";
            }
            else if(elementIsInBrackets)
            {
                bracketString += segment;
            }
            else if(!elementIsInBrackets)
            {
                parts.push_back(segment);
            }
        }

        if (parts.size() == 3) {
            // Handle cases where barcode might contain underscores
            std::string pattern = parts[0];
            std::string position = parts[1];
            std::string barcode  = parts[2];

            out << stripQuotes(pattern) << '\t' << position << '\t' << barcode << '\t' << "INSERTION" << '\t' << value << '\n';
        } else {
            std::cerr << "Invalid key format for statistics: " << key << '\n';
        }
    }
     //write deletions
     for (const auto& [key, value] : deletions) {
        std::stringstream ss(key);
        std::string segment;
        std::vector<std::string> parts;

        bool elementIsInBrackets = false;
        std::string bracketString = "";
        while (std::getline(ss, segment, '_')) 
        {
            if(!elementIsInBrackets && segment.at(0) == '\'')
            {
                if(segment.back() == '\'')
                {
                    parts.push_back(segment);
                }
                else
                {
                    elementIsInBrackets = true;
                    bracketString += segment + "_";
                }
            }
            else if(elementIsInBrackets && segment.back() == '\'')
            {
                assert(elementIsInBrackets);
                elementIsInBrackets = false;
                bracketString += segment;
                parts.push_back(bracketString);
                bracketString = "";
            }
            else if(elementIsInBrackets)
            {
                bracketString += segment;
            }
            else if(!elementIsInBrackets)
            {
                parts.push_back(segment);
            }
        }

        if (parts.size() == 3) {
            // Handle cases where barcode might contain underscores
            std::string pattern = parts[0];
            std::string position = parts[1];
            std::string barcode = parts[2];

            out << stripQuotes(pattern) << '\t' << position << '\t' << barcode << '\t' << "DELETION" << '\t' << value << '\n';
        } else {
            std::cerr << "Invalid key format: " << key << '\n';
        }
    }
     //write substitutions
     for (const auto& [key, value] : substitutions) {
        std::stringstream ss(key);
        std::string segment;
        std::vector<std::string> parts;

        bool elementIsInBrackets = false;
        std::string bracketString = "";
        while (std::getline(ss, segment, '_')) 
        {
            if(!elementIsInBrackets && segment.at(0) == '\'')
            {
                if(segment.back() == '\'')
                {
                    parts.push_back(segment);
                }
                else
                {
                    elementIsInBrackets = true;
                    bracketString += segment + "_";
                }
            }
            else if(elementIsInBrackets && segment.back() == '\'')
            {
                assert(elementIsInBrackets);
                elementIsInBrackets = false;
                bracketString += segment;
                parts.push_back(bracketString);
                bracketString = "";
            }
            else if(elementIsInBrackets)
            {
                bracketString += segment;
            }
            else if(!elementIsInBrackets)
            {
                parts.push_back(segment);
            }
        }

        if (parts.size() == 3) {
            // Handle cases where barcode might contain underscores
            std::string pattern = parts[0];
            std::string position = parts[1];
            std::string barcode = parts[2];

            out << stripQuotes(pattern) << '\t' << position << '\t' << barcode << '\t' << "SUBSTITUTION" << '\t' << value << '\n';
        } else {
            std::cerr << "Invalid key format: " << key << '\n';
        }
    }

    out.close();
}

void DemultiplexingStats::write_mm_number(const std::string& outputFile)
{
    //Find the maximum vector size (max number of mismatches)
    size_t maxMismatches = 0;
    for (const auto& [key, vec] : mismatchNumber) 
    {
        if (vec.size() > maxMismatches) {
            maxMismatches = vec.size();
        }
    }

    std::ofstream out(outputFile);
    if (!out) {
        std::cerr << "Could not open barcode-mismatche file for writing: " << outputFile << "\n";
        return;
    }

    //Write header
    out << "PATTERN\tPOSITION\tBARCODE";
    for (size_t i = 0; i < maxMismatches; ++i) {
        out << "\t" << i << "MM";
    }
    out << "\n";

    //Write each line
    for (const auto& [key, vec] : mismatchNumber) 
    {
        std::stringstream ss(key);
        std::string segment;
        std::vector<std::string> parts;

        bool elementIsInBrackets = false;
        std::string bracketString = "";
        while (std::getline(ss, segment, '_')) 
        {
            if(!elementIsInBrackets && segment.at(0) == '\'')
            {
                if(segment.back() == '\'')
                {
                    parts.push_back(segment);
                }
                else
                {
                    elementIsInBrackets = true;
                    bracketString += segment + "_";
                }
            }
            else if(elementIsInBrackets && segment.back() == '\'')
            {
                assert(elementIsInBrackets);
                elementIsInBrackets = false;
                bracketString += segment;
                parts.push_back(bracketString);
                bracketString = "";
            }
            else if(elementIsInBrackets)
            {
                bracketString += segment;
            }
            else if(!elementIsInBrackets)
            {
                parts.push_back(segment);
            }
        }

        if (parts.size() == 3) {
            std::string pattern = parts[0];
            std::string position = parts[1];
            std::string barcode = parts[2];

            out << stripQuotes(pattern) << '\t' << position << '\t' << barcode;

            // Now write mismatch counts
            for (size_t i = 0; i < maxMismatches; ++i) {
                if (i < vec.size()) {
                    out << '\t' << vec[i];
                } else {
                    out << "\t-";  // Fill missing values with -
                }
            }
            out << "\n";
        } else {
            std::cerr << "Invalid key format: (" << key << ") in barcode-mismatche file for writing: " << outputFile << "\n";
        }
    }

    out.close();
}


void DemultiplexingStats::write_last_mapped_position(const std::string& outputFile) 
{
    std::ofstream out(outputFile);
    if (!out) {
        std::cerr << "Could not open file for writing: " << outputFile << "\n";
        return;
    }

    // Header
    out << "PATTERN\tPOSITION\tREAD_DIRECTION\tCOUNT\n";

    //WRITE FORWARD
    for (const auto& [key, value] : failedLinesMappingFw) {
        std::stringstream ss(key);
        std::string pattern, position;

        // Split into pattern and position
        std::string segment;
        std::vector<std::string> parts;

        bool elementIsInBrackets = false;
        std::string bracketString = "";
        while (std::getline(ss, segment, '_')) 
        {
            if(!elementIsInBrackets && segment.at(0) == '\'')
            {
                if(segment.back() == '\'')
                {
                    parts.push_back(segment);
                }
                else
                {
                    elementIsInBrackets = true;
                    bracketString += segment + "_";
                }
            }
            else if(elementIsInBrackets && segment.back() == '\'')
            {
                assert(elementIsInBrackets);
                elementIsInBrackets = false;
                bracketString += segment;
                parts.push_back(bracketString);
                bracketString = "";
            }
            else if(elementIsInBrackets)
            {
                bracketString += segment;
            }
            else if(!elementIsInBrackets)
            {
                parts.push_back(segment);
            }
        }

        if(parts.size() == 2)
        {
            pattern = parts[0];
            position = parts[1];
            out << stripQuotes(pattern) << '\t' << position << '\t' << "FORWARD_READ" << '\t' << value << '\n';
        }
        else
        {
            std::cerr << "Invalid key: " << key << " in statistics for last mapped barcode position\n";
        }

    }

    //WRITE REVERSE
    for (const auto& [key, value] : failedLinesMappingRv) {
        std::stringstream ss(key);
        std::string pattern, position;

        // Split into pattern and position
        std::string segment;
        std::vector<std::string> parts;

        bool elementIsInBrackets = false;
        std::string bracketString = "";
        while (std::getline(ss, segment, '_')) 
        {
            if(!elementIsInBrackets && segment.at(0) == '\'')
            {
                if(segment.back() == '\'')
                {
                    parts.push_back(segment);
                }
                else
                {
                    elementIsInBrackets = true;
                    bracketString += segment + "_";
                }
            }
            else if(elementIsInBrackets && segment.back() == '\'')
            {
                assert(elementIsInBrackets);
                elementIsInBrackets = false;
                bracketString += segment;
                parts.push_back(bracketString);
                bracketString = "";
            }
            else if(elementIsInBrackets)
            {
                bracketString += segment;
            }
            else if(!elementIsInBrackets)
            {
                parts.push_back(segment);
            }
        }

        if(parts.size() == 2)
        {
            pattern = parts[0];
            position = parts[1];
            out << stripQuotes(pattern) << '\t' << position << '\t' << "REVERSE_READ" << '\t' << value << '\n';
        }
        else
        {
            std::cerr << "Invalid key: " << key << " in statistics for last mapped barcode position\n";
        }

    }

    out.close();
}

void DemultiplexingStats::write(const std::string& directory, const std::string& prefix, const int patternNumber)
{
    //create file names (potentially with prefix)
    std::string barcodeMismatchNumberFile = "Quality_numberMM.txt";
    std::string barcodeMismatchTypeFile = "Quality_typeMM.txt";
    std::string barcodeLastPosMappedFile = "Quality_lastPositionMapped.txt";

    if(prefix != "")
    {
        barcodeMismatchNumberFile = prefix + '_' + barcodeMismatchNumberFile;
        barcodeMismatchTypeFile = prefix + '_' + barcodeMismatchTypeFile;
        barcodeLastPosMappedFile = prefix + '_' + barcodeLastPosMappedFile;
    }

    //remove files if they already exist
    std::string barcodeMismatchNumber = directory + "/" + barcodeMismatchNumberFile;
    std::string barcodeMismatchType = directory + "/" + barcodeMismatchTypeFile;
    std::string barcodeLastPosMapped = directory + "/" + barcodeLastPosMappedFile;

    // remove outputfile if it exists
    std::remove(barcodeMismatchNumber.c_str());
    std::remove(barcodeMismatchType.c_str());
    std::remove(barcodeLastPosMapped.c_str());

    //fill the number of MM
    write_mm_number(barcodeMismatchNumber);

    //fill the type of MM
    write_mm_types(barcodeMismatchType);

    //fill last mapped positions
    if(patternNumber == 1)
    {
        write_last_mapped_position(barcodeLastPosMapped);
    }
}