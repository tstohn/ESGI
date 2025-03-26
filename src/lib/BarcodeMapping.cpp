
#include "BarcodeMapping.hpp"

bool check_if_seq_too_short(const int& offset, const std::string& seq)
{
    if(offset >= seq.length())
    {
        return true;
    }

    return false;
}

//initialize the dictionary of mismatches for each possible barcode
// the vector in the dict has length "mismatches + 2" for each barcode
// one entry for zero mismatches, eveery number from, 1 to mismatches and 
//one for more mismatches than the max allowed number of mismathces
template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::initializeStats()
{
    for(BarcodePatternPtr patternPtr : *barcodePatternList)
    {
        for(BarcodeVector::iterator patternItr = patternPtr->begin(); 
            patternItr < patternPtr->end(); 
            ++patternItr)
        {
            int mismatches = (*patternItr)->mismatches;
            const std::vector<std::string> patterns = (*patternItr)->get_patterns();
            for(const std::string& pattern : patterns)
            {
                std::vector<int> mismatchVector(mismatches + 2, 0);
                stats.mapping_dict.insert(std::make_pair(pattern, mismatchVector));
            }
        }   
    }
}

//parse a new file with barcodes and writes all barcodes into a vector
template <typename MappingPolicy, typename FilePolicy>
std::vector<std::string> Mapping<MappingPolicy, FilePolicy>::parse_variable_barcode_file(const std::string& barcodeFile)
{
    std::vector<std::string> barcodes;
    try
    {
        std::ifstream barcodeFileStream;
        barcodeFileStream.open(barcodeFile);
        std::stringstream strStream;
        strStream << barcodeFileStream.rdbuf();
        std::string barcodeListString = strStream.str();

        std::string delimiter = ",";
        int pos = 0;
        while ((pos = barcodeListString.find(delimiter)) != std::string::npos) {
            std::string seq = barcodeListString.substr(0, pos);
            barcodeListString.erase(0, pos + 1);
            for (char const &c: seq) {
                if(!(c=='A' || c=='T' || c=='G' || c=='C' ||
                        c=='a' || c=='t' || c=='g' || c=='c'))
                        {
                            std::cout << seq << "\n";
   
                        std::cerr << "PARAMETER ERROR: a barcode sequence in " << barcodeFile << " is not a base (A,T,G,C)\n";
                        if(c==' ' || c=='\t' || c=='\n')
                        {
                            std::cerr << "PARAMETER ERROR: Detected a whitespace in following sequence[" << seq << "], pls. remove it to continue!\n";
                            if(c=='\n')
                            {
                                std::cerr << "Looks like you forgot to remove a newline at the end of the file?\n";
                            }
                        }
                        exit(1);
                        }
            }
            barcodes.push_back(seq);
        }
        for (char const &c: barcodeListString) {
            if(!(c=='A' || c=='T' || c=='G' || c=='C' ||
                    c=='a' || c=='t' || c=='g' || c=='c'))
                    {
                    std::cerr << "PARAMETER ERROR: a barcode sequence in " << barcodeFile << " is not a base (A,T,G,C)\n";
                    if(c==' ' || c=='\t' || c=='\n')
                    {
                        std::cerr << "PARAMETER ERROR: Detected a whitespace in following sequence[" << barcodeListString << "], pls. remove it to continue!\n";
                        if(c=='\n')
                        {
                            std::cerr << "Looks like you forgot to remove a newline at the end of the file?\n";
                        }
                    }
                    exit(1);
                    }
        }
        barcodes.push_back(barcodeListString);
            
        barcodeFileStream.close();
    }
    catch(std::exception& e)
    {
        std::cerr << "Detected a file name for variable barcode patterns but file could not be parsed! File: " << barcodeFile << "\n";
        std::cerr << e.what() << std::endl;
        exit(1);
    }

    return(barcodes);
}

//parse a single line of the pattern file like: [GACTTCAG][15X][barcodes.txt][DNA]
template <typename MappingPolicy, typename FilePolicy>
std::pair<std::string, std::vector<std::string> > Mapping<MappingPolicy, FilePolicy>::parse_pattern_line(std::string& line, int number)
{
    std::vector<std::string> patterns;
    std::string name;
    try{
        //if present parse the name of the line
        size_t pos = line.find(":[");
        if (pos != std::string::npos) 
        {
            //make sure this patterns only exists once
            if (line.find(":[", pos + 1) != std::string::npos) 
            {
                std::cerr << "The substring <:[> appears more than once. This is invalid as <:> should only be used to" <<
                " name patterns like, e.g., <name:[GACGTA][pattern.txt]>... " << std::endl;
                exit(EXIT_FAILURE);     
            }

            //parse name
            name = line.substr(0, pos);
            line.erase(0, pos + 1);
        }
        else
        {
            name = "PATTERN_" + std::to_string(number);
        }

        // parse one pattern line: e.g.: [ACGTTCAG][15X][file1.txt]
        std::string delimiter = "]";
        const char delimiter2 = '[';
        pos = 0;
        std::string seq;
        //PARSE PATTERN
        while ((pos = line.find(delimiter)) != std::string::npos) {
            seq = line.substr(0, pos);
            line.erase(0, pos + 1);
            if(seq.at(0) != delimiter2)
            {
                std::cerr << "PARAMETER ERROR: Wrong barcode patter parameter, it looks like a \'[\' is missing\n";
                exit(1);
            }
            seq.erase(0, 1);
            
            patterns.push_back(seq);
        }
    }
    catch(std::exception& e)
    {
        std::cerr << "PARAMETER ERROR: Could not parse following line in the barcode pattern file: " << line << "\n";
        std::cerr << e.what() << std::endl;
        exit(1);
    }

    std::pair<std::string, std::vector<std::string>> parsedLinePair;

    // Initialize the pair
    parsedLinePair.first = name;
    parsedLinePair.second = patterns;

    return(parsedLinePair);
}

//parsing a new pattern line of the patterns file and store it in a vector
//it also parses names for lines if present:
//GUIDELINE:[GCAGCT][10X][GUIDES.txt]
template <typename MappingPolicy, typename FilePolicy>
std::vector<std::pair<std::string, std::vector<std::string>>> Mapping<MappingPolicy, FilePolicy>::parse_pattern_file(const std::string& patternFile)
{

    std::vector<std::pair<std::string, std::vector<std::string>>> patternList;

    //iterate through rows in pattern file: whenever we encounter a new file of barcodes parse it
    std::ifstream patternFileStream(patternFile);
    if(!patternFileStream.is_open())
    {
        std::cerr << "Opening File with barcode patterns failed. File name: " << patternFile;
        exit(EXIT_FAILURE);
    }
    std::string line;
    int lineNumber = 0;
    while(std::getline(patternFileStream, line))
    {
        patternList.emplace_back(parse_pattern_line(line, lineNumber));
        ++lineNumber;
    }

    return(patternList);
}

//parse a new lnie in mismatch file
template <typename MappingPolicy, typename FilePolicy>
std::vector<std::vector<int>> Mapping<MappingPolicy, FilePolicy>::parse_mismatch_file(const std::string& mismatchFile)
{
    std::vector<std::vector<int>> mismatchList;

    std::ifstream mismatchFileStream(mismatchFile);
    if(!mismatchFileStream.is_open())
    {
        std::cerr << "Opening File with mismatches per pattern failed. File name: " << mismatchFile;
        exit(EXIT_FAILURE);
    }
    std::string mismatchLine;

    //parse every line in the mismatch file, specific for the corresponding pattern line in patternFile
    std::string delimiter = ",";
    while(std::getline(mismatchFileStream, mismatchLine))
    {
        std::vector<int> mismatches;
        int pos = 0;
        while ((pos = mismatchLine.find(delimiter)) != std::string::npos) 
        {
            std::string seq = mismatchLine.substr(0, pos);
            mismatchLine.erase(0, pos + 1);
            mismatches.push_back(std::stoi(seq));
        }
        mismatches.push_back(std::stoi(mismatchLine));

        //push the vector of mismatches into the mismatchList which stores all vectors of mismatches
        //for every pattern line
        mismatchList.push_back(mismatches);
    }

    return(mismatchList);
}
        
//creates a barcodePattern(which is essentially a list of barcodes wrapped around a functional class)
template <typename MappingPolicy, typename FilePolicy>
BarcodePatternPtr Mapping<MappingPolicy, FilePolicy>::create_barcodeVector_from_patternLine(
    const std::vector<std::string>& barcodeList, 
    const std::vector<int>& mismatchList, 
    const std::string& patternName,
    std::unordered_map<std::string, std::vector<std::string>>& fileToBarcodesMap)
{

    BarcodeVector barcodeVector;

    bool containsDNA = false;
    //iterate through all barcodes in the pattern line
    //parse the pattern and immediately create the Barcode
    for(int barcodeIdx = 0; barcodeIdx < barcodeList.size(); ++barcodeIdx)
    {
        std::string patternElement = barcodeList.at(barcodeIdx);
        int barcodeLength = 0;
        bool barcodeFound = false;

        // constant barcode
        bool isConstant = true;
        for (char const &c: patternElement) 
        {
            if(c!='A' && c!='T' && c!='G' && c!='C' &&
            c!='a' && c!='t' && c!='g' && c!='c')
            {
                isConstant = false;
                break;
            }
        }
        if(isConstant)
        {
            ConstantBarcode barcode(patternElement, mismatchList.at(barcodeIdx));
            std::shared_ptr<ConstantBarcode> barcodePtr(std::make_shared<ConstantBarcode>(barcode));
            barcodeVector.push_back(barcodePtr);
            barcodeFound = true;
        }

        //variable barcode file
        //check if the string is in our file->barcodes map
        //if not check if we can parse it
        bool isVariable = false;
        //if the file was already parsed we know this must be a variable barcode
        if(fileToBarcodesMap.find(patternElement) != fileToBarcodesMap.end())
        {
            isVariable = true;
        }
        else //otherwise try to open it as a file
        {
            std::ifstream possibleFileStream(patternElement);
            if(possibleFileStream.good())
            {
                //if possible store it in the map of files -> barcode vector
                fileToBarcodesMap[patternElement] = parse_variable_barcode_file(patternElement);
                isVariable = true;
            }
        }
        //create Variable Barcode from this data
        if(isVariable)
        {
            VariableBarcode barcode(fileToBarcodesMap.at(patternElement), patternElement, mismatchList.at(barcodeIdx));
            std::shared_ptr<VariableBarcode> barcodePtr(std::make_shared<VariableBarcode>(barcode));
            barcodeVector.push_back(barcodePtr);
            barcodeFound = true;
        }

        //random element(e.g., UMI)
        //format: [15X], before we must check if the prefix-digit is not part of a file name (done above for variable barcode file)
        if(std::isdigit(patternElement[0]))
        {
            size_t pos = 0;
            barcodeLength = std::stoi(patternElement, &pos);
            std::string seq = patternElement.substr(pos);
            if(seq == "X")
            {
                //For random sequences we MUST specific the length of this sequence
                WildcardBarcode barcode(mismatchList.at(barcodeIdx), patternElement, barcodeLength);
                std::shared_ptr<WildcardBarcode> barcodePtr(std::make_shared<WildcardBarcode>(barcode));
                barcodeVector.push_back(barcodePtr);
                barcodeFound = true;
            }
            else
            {
                //we already checked above if it is an existing barcode file
                std::cerr << "Encountered an error in barode pattern: " << patternElement << ".\n" <<
                "If a pattern has a digit as a prefix it must ether be part of a barcode-file name " <<
                "or it must be a random sequence like UMIs in the format <NUMBER><X> with a number followed by a sinlge capital X." << 
                "For example [15X]\n";
                exit(1);
            }
        }

        //DNA element
        if(patternElement == "DNA")
        {
            DNABarcode barcode(mismatchList.at(barcodeIdx));
            std::shared_ptr<DNABarcode> barcodePtr(std::make_shared<DNABarcode>(barcode));
            barcodeVector.push_back(barcodePtr);
            barcodeFound = true;
            containsDNA = true;
        }

        //stop pattern: [*] (only map up to here)
        if(patternElement == "*")
        {
            StopBarcode barcode(patternElement, mismatchList.at(barcodeIdx));
            std::shared_ptr<StopBarcode> barcodePtr(std::make_shared<StopBarcode>(barcode));
            barcodeVector.push_back(barcodePtr);
            barcodeFound = true;
        }

        //seperator pattern: [-] (ether stop here, or if DNA extract all till the read ends)
        if(patternElement == "-")
        {
            ReadSeperatorBarcode barcode(patternElement, mismatchList.at(barcodeIdx));
            std::shared_ptr<ReadSeperatorBarcode> barcodePtr(std::make_shared<ReadSeperatorBarcode>(barcode));
            barcodeVector.push_back(barcodePtr);
            barcodeFound = true;
        }

        //if nothing can be created from it, return false and throw an error
        if(!barcodeFound)
        {
            std::cerr << "The pattern element [" << patternElement <<"] is not a valid barcode element." <<
            "Read the manual again to make sure it matches, e.g.: a constant barcode (only A,G,T,C), a file which contains all possible barcodes" <<
            "at this position, a random barcode like UMIs with <number><X> like <15X>, or a simple string such as <DNA> for a DNA/ RNA sequence," <<
            "or <*> for a stop position (to seperate FW, RV reads, or simply exclude a sequence in the middle).";
            exit(1);
        }

    }

    BarcodeVectorPtr barcodeVectorPtr = std::make_shared<BarcodeVector>(barcodeVector);
    BarcodePatternPtr patternPtr = std::make_shared<BarcodePattern>(BarcodePattern(containsDNA, patternName, barcodeVectorPtr));

    return(patternPtr);
}

//new fucntion to overwrite old function which was misconstructed
template <typename MappingPolicy, typename FilePolicy>
bool Mapping<MappingPolicy, FilePolicy>::generate_barcode_patterns(const input& input)
{
    
    //a map of file names to the actual barcodes in this file. This map is updated while itereating through the pattern file
    //whenever an unknown file is encountered
    std::unordered_map<std::string, std::vector<std::string>> fileToBarcodesMap;

    //only arsing, no quality check
    std::vector<std::pair<std::string, std::vector<std::string>>> patternList = parse_pattern_file(input.barcodePatternsFile);
    std::vector<std::vector<int>> mismatchList = parse_mismatch_file(input.mismatchFile);

    // assert that the number of mismatches and patterns are the same
    if(patternList.size() != mismatchList.size())
    {
        std::cerr << "Number of barcode pattern lines and mismatch lines does not match. Please correct in parameters.\n";
        exit(1);
    }

    for(int i = 0; i < patternList.size(); i++)
    {
        if(patternList.at(i).second.size() != mismatchList.at(i).size())
        {
            std::cerr << "Number of barcode patterns and mismatches does not match for line " << std::to_string(i) << ". Please correct in parameters.\n";
            exit(EXIT_FAILURE);
        }
    }
    //make sure pattern names are unique
    std::unordered_map<std::string, int> countMap;
    for(int i = 0; i < patternList.size(); i++)
    {
        countMap[patternList.at(i).first]++;
    }
    for (const auto& [key, count] : countMap) 
    {
        if (count > 1) {
            std::cout << "Duplicate pattern names in pattern file (-p), e.g., " << key << " \n";
            std::cout << "Remove duplicate names, or remove names completely (a number as name will be given in order of the patterns) \n";           
            exit(EXIT_FAILURE);
        }
    }

    // parse all files with variables barcodes (guides, BC1, BC2, ..., ABs)
    // parse mismatches, pattern, lengths of patterns (e.g. UMI legnth 15)
    for(int i = 0; i < patternList.size(); i++)
    {
        barcodePatternList->emplace_back(create_barcodeVector_from_patternLine(patternList.at(i).second, mismatchList.at(i), patternList.at(i).first, fileToBarcodesMap));
    }

    return true;
}

//apply all BarcodePatterns to a single line to generate a vector strings (the Barcodemapping)
//realBarcodeMap contains the actual string in the seauence that gets a barcode assigned
bool MapEachBarcodeSequentiallyPolicy::split_line_into_barcode_patterns(const std::pair<fastqLine, fastqLine>& seq, 
                                        DemultiplexedLine& demultiplexedLine, const input& input, 
                                        BarcodePatternPtr barcodePatterns,
                                        fastqStats& stats)
{

    //iterate over BarcodeMappingVector
    int offset = 0;
    int score_sum = 0;

    int old_offset = offset;
    unsigned int wildCardToFill = 0;
    int wildCardLength = 0, differenceInBarcodeLength = 0;
    for(BarcodeVector::iterator patternItr = barcodePatterns->begin(); 
        patternItr < barcodePatterns->end(); 
        ++patternItr)
    {

        //if we have a wildcard skip this matching, we match again the next sequence
        if((*patternItr)->is_wildcard())
        {
            if(!wildCardToFill)
            {
                old_offset = offset;
            }
            wildCardLength = (*patternItr)->length;
            offset += wildCardLength;
            wildCardToFill += 1;
            continue;
        }
        else if((*patternItr)->is_stop())
        {
            //stop here: we do not continue mapping after stop barcode [*]
            break;
        }
        else if((*patternItr)->is_dna())
        {
            //make sure DNA is the last part of a sequence
            // (we do not know length of DNA and can not map barcodes on both sides of the DNA)
            if(!(*(patternItr+1))->is_read_end())
            {
                std::cerr << "After a DNA pattern a read-end pattern [-] must follow, since we do not know the length of the DNA pattern.\n\
                 Therefore a valid pattern would be ether ... [DNA][-]... or ...[-][DNA]...\n\
                Please adjust pattern file accordingly.";
                exit(EXIT_FAILURE);
            }

            //set dna, quality (name was already set when getting next line to the name of forward read ONLY)
            demultiplexedLine.dna = seq.first.line.substr(offset, seq.first.line.length()); //extract DNA fragment
            demultiplexedLine.dnaQuality = seq.first.quality.substr(offset, seq.first.line.length()); //extract DNA fragment
            demultiplexedLine.containsDNA = true;
        }
        else if((*patternItr)->is_read_end())
        {
            //stop here: we do not continue mapping when the read ends
            break;
        }

        //for every barcodeMapping element find a match
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
        int start=0, end=0, score = 0;
        bool startCorrection = false;
        if(wildCardToFill){startCorrection = true;} // sart correction checks if we have to move our mapping window to the 5' direction
        // could happen in the case of deletions in the UMI sequence...
        bool seqToShort = check_if_seq_too_short(offset, seq.first.line);
        if(seqToShort)
        {
            ++stats.noMatches;
            return false;
        }

        if(!(*patternItr)->match_pattern(seq.first.line, offset, start, end, score, barcode, differenceInBarcodeLength, startCorrection, false))
        {
            ++stats.noMatches;
            return false;
        }
        
        std::string sequence = seq.first.line.substr(offset + start, end-start);
        offset += end;
        score_sum += score;

        assert(barcode != "");
        if(input.writeStats)
        {
            //add barcode data to statistics dictionary
            std::lock_guard<std::mutex> guard(*stats.statsLock);
            int dictvectorIndex = ( (score <= (*patternItr)->mismatches) ? (score) : ( ((*patternItr)->mismatches) + 1) );
            ++stats.mapping_dict[barcode].at(dictvectorIndex);
        }
        
        //squeeze in the last wildcard match if there was one 
        //(we needed both matches of neighboring barcodes to define wildcard boundaries)
        if(wildCardToFill)
        {
            startCorrection = false;
            std::string oldWildcardMappedBarcode = seq.first.line.substr(old_offset, (offset+start-end) - old_offset);
            //barcodeMap.emplace_back(uniqueChars.  (oldWildcardMappedBarcode));
            unsigned int lastWildcardEnd = 0; //in case we have several wildcards and need to know old offset
            // to subset new string
            while(wildCardToFill != 0)
            {
                BarcodeVector::iterator wildCardIt = patternItr;
                int pos = -1 * wildCardToFill;
                std::advance(wildCardIt, pos);
                int currentWildcardLength = (*wildCardIt)->length;
                //in case of deletions in UMI, we remove nucleotides from the last wildcard sequence
                //for the last UMI sequence, we take the whole sequence (in case of insertions)

                if( (lastWildcardEnd + currentWildcardLength) > oldWildcardMappedBarcode.length() || wildCardToFill == 1) 
                {
                    currentWildcardLength = oldWildcardMappedBarcode.length() - lastWildcardEnd;
                }

                std::string wildCardString = oldWildcardMappedBarcode.substr(lastWildcardEnd, currentWildcardLength);

                demultiplexedLine.barcodeList.push_back(wildCardString);
                wildCardToFill -= 1;

                lastWildcardEnd += currentWildcardLength; // add the lengths of barcodes to the next offset position
            }
        }
        //add this match to the BarcodeMapping
        //barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        demultiplexedLine.barcodeList.push_back(barcode);
    }
    //if the last barcode was a WIldcard that still has to be added:
    //BE CAREFUL: for now this means there can be ONLY ONE UMI at the end of a sequence
    if(wildCardToFill)
    {
        assert(wildCardToFill <= 1);
        wildCardToFill = 0; // unnecessary, still left to explicitely set to false
        std::string oldWildcardMappedBarcode = seq.first.line.substr(old_offset, seq.first.line.length() - old_offset);
        //barcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
        demultiplexedLine.barcodeList.push_back(oldWildcardMappedBarcode);
    }

    if(score_sum == 0)
    {
        ++stats.perfectMatches;
    }
    else
    {
        ++stats.moderateMatches;
    }

    return true;
}

bool MapEachBarcodeSequentiallyPolicyPairwise::map_forward(const fastqLine& seq, const input& input, 
                                                           BarcodePatternPtr barcodePatterns,
                                                           fastqStats& stats,
                                                           DemultiplexedLine& demultiplexedLine,
                                                           uint& barcodePosition,
                                                           int& score_sum)
{

    //iterate over BarcodeMappingVector
    int offset = 0;

    int old_offset = offset;
    unsigned int wildCardToFill = 0;
    int wildCardLength = 0, differenceInBarcodeLength = 0;
    for(BarcodeVector::iterator patternItr = barcodePatterns->begin(); 
        patternItr < barcodePatterns->end(); 
        ++patternItr)
    {
        //if we have a wildcard skip this matching, we match again the next sequence
        if((*patternItr)->is_wildcard())
        {
            if(!wildCardToFill)
            {
                old_offset = offset;
            }
            wildCardLength = (*patternItr)->length;
            offset += wildCardLength;
            wildCardToFill += 1;
            continue;
        }
        else if((*patternItr)->is_stop())
        {
            //the stop positions is also a found position (count only for fw)
            ++barcodePosition; //increase the count of found positions
            //stop here: we do not continue mapping after stop barcode [*]
            return true;
        }
        else if((*patternItr)->is_dna())
        {
            //make sure DNA is the last part of a sequence
            // (we do not know length of DNA and can not map barcodes on both sides of the DNA)
            if(!(*(patternItr+1))->is_read_end())
            {
                std::cerr << "After a DNA pattern a read-end pattern [-] must follow, since we do not know the length of the DNA pattern.\n\
                 Therefore a valid pattern would be ether ... [DNA][-]... or ...[-][DNA]...\n\
                Please adjust pattern file accordingly.";
                exit(EXIT_FAILURE);
            }
            ++barcodePosition; //increase the count of found positions for DNA

            //set dna, quality (name was already set when getting next line to the name of forward read ONLY)
            demultiplexedLine.dna = seq.line.substr(offset, seq.line.length()); //extract DNA fragment
            demultiplexedLine.dnaQuality = seq.quality.substr(offset, seq.line.length()); //extract DNA fragment
            demultiplexedLine.containsDNA = true;
            continue;
        }
        else if((*patternItr)->is_read_end())
        {
            ++barcodePosition; //read-end is a valid found barcode
            //stop here: we do not continue mapping when the read ends
            return true;
        }

        //for every barcodeMapping element find a match
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
        int start=0, end=0, score = 0;
        bool startCorrection = false;
        if(wildCardToFill){startCorrection = true;} // sart correction checks if we have to move our mapping window to the 5' direction
        
        // in case the sequence ends in the middle of the next barcode
        bool seqToShort = check_if_seq_too_short(offset, seq.line);
        if(seqToShort)
        {
            return true;
        }

        //if we did not match a pattern
        if(!(*patternItr)->match_pattern(seq.line, offset, start, end, score, barcode, differenceInBarcodeLength, startCorrection, false))
        {
            return false;
        }

        //set the length difference after barcode mapping
        //for the case of insertions inside the barcode sequence set the difference explicitely to zero 
        //(we only focus on deletions that we can not distinguish from substitutions)
        //differenceInBarcodeLength = barcode.length() - (end);
        if(differenceInBarcodeLength<0){differenceInBarcodeLength=0;}
        std::string sequence = seq.line.substr(offset + start, end-start);
        offset += end;
        score_sum += score;

        assert(barcode != "");
        if(input.writeStats)
        {
            //add barcode data to statistics dictionary
            std::lock_guard<std::mutex> guard(*stats.statsLock);
            int dictvectorIndex = ( (score <= (*patternItr)->mismatches) ? (score) : ( ((*patternItr)->mismatches) + 1) );
            ++stats.mapping_dict[barcode].at(dictvectorIndex);
        }
        
        //squeeze in the last wildcard match if there was one 
        //(we needed both matches of neighboring barcodes to define wildcard boundaries)
        if(wildCardToFill)
        {
            startCorrection = false;
            std::string oldWildcardMappedBarcode = seq.line.substr(old_offset, (offset+start-end) - old_offset);
            //barcodeMap.emplace_back(uniqueChars.  (oldWildcardMappedBarcode));
            unsigned int lastWildcardEnd = 0; //in case we have several wildcards and need to know old offset
            // to subset new string
            while(wildCardToFill != 0)
            {
                BarcodeVector::iterator wildCardIt = patternItr;
                int pos = -1 * wildCardToFill;
                std::advance(wildCardIt, pos);
                int currentWildcardLength = (*wildCardIt)->length;
                //in case of deletions in UMI, we remove nucleotides from the last wildcard sequence
                //for the last UMI sequence, we take the whole sequence (in case of insertions)
                if( (lastWildcardEnd + currentWildcardLength) > oldWildcardMappedBarcode.length() || wildCardToFill == 1) 
                {
                    currentWildcardLength = oldWildcardMappedBarcode.length() - lastWildcardEnd;
                }

                std::string wildCardString = oldWildcardMappedBarcode.substr(lastWildcardEnd, currentWildcardLength);
                demultiplexedLine.barcodeList.push_back(wildCardString);
                wildCardToFill -= 1;

                lastWildcardEnd += currentWildcardLength; // add the lengths of barcodes to the next offset position
                ++barcodePosition; //increase the count of found positions
            }
        }

        //add this match to the BarcodeMapping
        //barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        demultiplexedLine.barcodeList.push_back(barcode);
        ++barcodePosition; //increase the count of found positions
    }
    //if the last barcode was a WIldcard that still has to be added
    if(wildCardToFill)
    {
        wildCardToFill = 0; // unnecessary, still left to explicitely set to false
        std::string oldWildcardMappedBarcode = seq.line.substr(old_offset, seq.line.length() - old_offset);
        //barcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
        demultiplexedLine.barcodeList.push_back(oldWildcardMappedBarcode);
        ++barcodePosition; //increase the count of found positions
    }

    return true;
}

bool MapEachBarcodeSequentiallyPolicyPairwise::map_reverse(const fastqLine& seq, const input& input, 
                                                           BarcodePatternPtr barcodePatterns,
                                                           fastqStats& stats,
                                                           DemultiplexedLine& demultiplexedLine,
                                                           uint& barcodePosition,
                                                           int& score_sum)
{
    //iterate over BarcodeMappingVector
    int offset = 0;

    int old_offset = offset;
    unsigned int wildCardToFill = 0;
    int wildCardLength = 0, differenceInBarcodeLength = 0;
    //iterate reverse through patterns
    for(BarcodeVector::reverse_iterator patternItr = barcodePatterns->rbegin(); 
        patternItr < barcodePatterns->rend(); 
        ++patternItr)
    {
        //if we have a wildcard skip this matching, we match again the next sequence
        if((*patternItr)->is_wildcard())
        {
            if(!wildCardToFill)
            {
                old_offset = offset;
            }
            wildCardLength = (*patternItr)->length;
            offset += wildCardLength;
            wildCardToFill += 1;
            continue;
        }
        else if((*patternItr)->is_stop())
        {
            //stop here: we do not continue mapping after stop barcode [*]
            return true;
        }
        else if((*patternItr)->is_dna())
        {
            //make sure DNA is the last part of a sequence
            // (we do not know length of DNA and can not map barcodes on both sides of the DNA)
            if(!(*(patternItr+1))->is_read_end())
            {
                std::cerr << "After a DNA pattern a read-end pattern [-] must follow, since we do not know the length of the DNA pattern.\n\
                 Therefore a valid pattern would be ether ... [DNA][-]... or ...[-][DNA]...\n\
                Please adjust pattern file accordingly.";
                exit(EXIT_FAILURE);
            }

            //set dna, quality (name was already set when getting next line to the name of forward read ONLY)
            demultiplexedLine.dna = seq.line.substr(offset, seq.line.length()); //extract DNA fragment
            demultiplexedLine.dnaQuality = seq.quality.substr(offset, seq.line.length()); //extract DNA fragment
            demultiplexedLine.containsDNA = true;
        }
        else if((*patternItr)->is_read_end())
        {
            //stop here: we do not continue mapping when the read ends
            break;
        }

        //for every barcodeMapping element find a match
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
        int start=0, end=0, score = 0;
        bool startCorrection = false;
        if(wildCardToFill){startCorrection = true;} // sart correction checks if we have to move our mapping window to the 5' direction
        
        // in case the sequence ends in the middle of the next barcode
        bool seqToShort = check_if_seq_too_short(offset, seq.line);
        if(seqToShort)
        {
            return true;
        }

        //map each pattern with reverse complement
        if(!(*patternItr)->match_pattern(seq.line, offset, start, end, score, barcode, differenceInBarcodeLength, startCorrection, true))
        {
            return false;
        }

        //set the length difference after barcode mapping
        //for the case of insertions inside the barcode sequence set the difference explicitely to zero 
        //(we only focus on deletions that we can not distinguish from substitutions)
        //we only substract the end, bcs we only consider barcode length differences at the end of the barcode
        //differenceInBarcodeLength = barcode.length() - (end);
        if(differenceInBarcodeLength<0){differenceInBarcodeLength=0;}
        std::string sequence = seq.line.substr(offset + start, end-start);
        offset += end;
        score_sum += score;

        assert(barcode != "");
        if(input.writeStats)
        {
            //add barcode data to statistics dictionary
            std::lock_guard<std::mutex> guard(*stats.statsLock);
            int dictvectorIndex = ( (score <= (*patternItr)->mismatches) ? (score) : ( ((*patternItr)->mismatches) + 1) );
            ++stats.mapping_dict[barcode].at(dictvectorIndex);
        }
        
        //squeeze in the last wildcard match if there was one 
        //(we needed both matches of neighboring barcodes to define wildcard boundaries)
        if(wildCardToFill)
        {
            startCorrection = false;
            std::string oldWildcardMappedBarcode = seq.line.substr(old_offset, (offset+start-end) - old_offset);

            //barcodeMap.emplace_back(uniqueChars.  (oldWildcardMappedBarcode));
            unsigned int lastWildcardEnd = 0; //in case we have several wildcards and need to know old offset
            // to subset new string
            while(wildCardToFill != 0)
            {
                BarcodeVector::reverse_iterator wildCardIt = patternItr;
                int pos = -1 * wildCardToFill;
                std::advance(wildCardIt, pos);
                int currentWildcardLength = (*wildCardIt)->length;
                //in case of deletions in UMI, we remove nucleotides from the last wildcard sequence
                //for the last UMI sequence, we take the whole sequence (in case of insertions)
                if( (lastWildcardEnd + currentWildcardLength) > oldWildcardMappedBarcode.length() || wildCardToFill == 1) 
                {
                    currentWildcardLength = oldWildcardMappedBarcode.length() - lastWildcardEnd;
                }

                std::string wildCardString = oldWildcardMappedBarcode.substr(lastWildcardEnd, currentWildcardLength);

                //in reverse mapping we have to make reverse complement of wildcard sequence
                std::string reverseComplimentbarcode = Barcode::generate_reverse_complement(wildCardString);

                demultiplexedLine.barcodeList.push_back(reverseComplimentbarcode);
                wildCardToFill -= 1;

                lastWildcardEnd += currentWildcardLength; // add the lengths of barcodes to the next offset position
                ++barcodePosition; //increase the count of found positions
            }
        }



        //add this match to the BarcodeMapping
        //barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        demultiplexedLine.barcodeList.push_back(barcode);
        ++barcodePosition; //increase the count of found positions
    }
    //if the last barcode was a WIldcard that still has to be added
    if(wildCardToFill)
    {
        wildCardToFill = 0; // unnecessary, still left to explicitely set to false
        std::string oldWildcardMappedBarcode = seq.line.substr(old_offset, seq.line.length() - old_offset);
        //barcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
        demultiplexedLine.barcodeList.push_back(oldWildcardMappedBarcode);
        ++barcodePosition; //increase the count of found positions
    }

    return true;
}

//demultiplexedLineFw is actually already the final barcodeList (stored in Demultiplexline), to not copy vectors of strings around but add in place
//if reverse read was DNA and contains that information we need to copy it into demultiplexedLineFw
bool MapEachBarcodeSequentiallyPolicyPairwise::combine_mapping(const BarcodePatternPtr& barcodePatterns,
                                                               DemultiplexedLine& demultiplexedLineFw,
                                                               const uint& barcodePositionFw,
                                                               const DemultiplexedLine& demultiplexedLineRv,
                                                               const uint& barcodePositionRv,
                                                               fastqStats& stats,
                                                               int& score_sum)
{
    //check that we span the whole sequence (except constant regions)
    int patternNum = barcodePatterns->size();

    //if the reverse read contains DNA co-y the information into the Fw demultiplexed line (the fw line is used as final class)
    if(demultiplexedLineRv.containsDNA)
    {
        demultiplexedLineFw.dna = demultiplexedLineRv.dna; //extract DNA fragment
        demultiplexedLineFw.dnaQuality = demultiplexedLineRv.dnaQuality; //extract DNA fragment
        demultiplexedLineFw.containsDNA = true;
    }

    //if positions are next to each other just return
    //in case of a stop pattern [*], we added +1 to the barcodePositionFw, so that barcodePositionFw+barcodePositionRv should be euqual to patternNum
    if(patternNum == (barcodePositionFw + barcodePositionRv))
    {
        for(std::vector<std::string>::const_reverse_iterator rvBarcodeIt = demultiplexedLineRv.barcodeList.crbegin(); rvBarcodeIt != demultiplexedLineRv.barcodeList.crend(); ++rvBarcodeIt)
        {
            demultiplexedLineFw.barcodeList.push_back({*rvBarcodeIt});
        }
    }
    //if they r overlapping, check all overlapping ones are the same
    if(patternNum < (barcodePositionFw + barcodePositionRv))
    {
        //assert all overlapping barcodes r the same
        int start = (patternNum - barcodePositionRv); // this is the idx of the first shared barcode
        int end = barcodePositionFw - 1; //this is the idnex of the last shared barcode
        int j = demultiplexedLineRv.barcodeList.size() - 1;
        for(int i = start; i <= end; ++i)
        {
            if(demultiplexedLineFw.barcodeList.at(i) != demultiplexedLineRv.barcodeList.at(j))
            {
                ++stats.noMatches;
                return false;
            }
            --j;
        }

        int overlap = (end - start) + 1;
        int missingPatternStartIdx = (demultiplexedLineRv.barcodeList.size() - overlap) -1;
        //combine them to one vector
        for(int i = missingPatternStartIdx; i >= 0; --i)
        {
            demultiplexedLineFw.barcodeList.push_back(demultiplexedLineRv.barcodeList.at(i));
        }
    }
    //if there r barcodePatterns missing in the middle, check they r only constant and add any
    //non-found-constant barcode
    if(patternNum > (barcodePositionFw + barcodePositionRv))
    {
        //check if missing barcodes are only constant
        int start = barcodePositionFw; // index of first missing barcode
        int end = (patternNum - barcodePositionRv) - 1; //last idnex that can be missing

        for(int i = start; i <= end; ++i)
        {
            if(!(barcodePatterns->barcodePattern->at(i)->is_constant()))
            {
                ++stats.noMatches;
                return false;
            }
            demultiplexedLineFw.barcodeList.push_back(barcodePatterns->barcodePattern->at(i)->get_patterns().at(0));
        }

        //add reverse patterns
        for(std::vector<std::string>::const_reverse_iterator rvBarcodeIt = demultiplexedLineRv.barcodeList.crbegin(); rvBarcodeIt != demultiplexedLineRv.barcodeList.crend(); ++rvBarcodeIt)
        {
            demultiplexedLineFw.barcodeList.push_back({*rvBarcodeIt});
        }
    }

    if(score_sum == 0)
    {
        ++stats.perfectMatches;
    }
    else
    {
        ++stats.moderateMatches;
    }

    return true;
}

bool MapEachBarcodeSequentiallyPolicyPairwise::split_line_into_barcode_patterns(const std::pair<fastqLine, fastqLine>& seq, 
                                        DemultiplexedLine& demultiplexedLine,
                                        const input& input,
                                        BarcodePatternPtr barcodePatterns,
                                        fastqStats& stats)
{

    int score_sum = 0;
    //map forward reads barcodeList contains the stored barcodes, 
    //barcodePosition is the psoiton of the last mapped barcode
    uint barcodePositionFw = 0;
    //for forward read we immediately add barcodes to demultiplexedLine.barcodeList, which is then extended in combine pattern, IF we find all patterns
    bool fwBool = map_forward(seq.first, input, barcodePatterns, stats, demultiplexedLine, barcodePositionFw, score_sum);

    DemultiplexedLine demultiplexedLineRv;
    uint barcodePositionRv = 0;

    bool rvBool = map_reverse(seq.second, input, barcodePatterns, stats, demultiplexedLineRv,barcodePositionRv, score_sum);

    //passing demultiplexedLine.barcodeList as the barcodeList in forward pattern
    bool pairwiseMappingSuccess = combine_mapping(barcodePatterns, demultiplexedLine, barcodePositionFw, demultiplexedLineRv, barcodePositionRv, stats, score_sum);

    return (pairwiseMappingSuccess);
}

template <typename MappingPolicy, typename FilePolicy>
bool Mapping<MappingPolicy, FilePolicy>::demultiplex_read(const std::pair<fastqLine, fastqLine>& seq, 
                                                          DemultiplexedLine& demultiplexedLine,
                                                          BarcodePatternPtr pattern,
                                                          const input& input, 
                                                          const unsigned long long& count, const unsigned long long& totalReadCount,
                                                          bool guideMapping)
{
    //split line into patterns (barcodeMap, barcodePatters, stats are passed as reference or ptr)
    //and can be read by each thread, "addValue" method for barcodeMap is thread safe also for concurrent writing
    bool result;

    //demultipelxed barcodes are stored in barcodeMap
    result = this->split_line_into_barcode_patterns(seq, demultiplexedLine, input, pattern, stats);

    //update status bar
    if(count%1000==0 && totalReadCount!=ULLONG_MAX && count<totalReadCount) //update at every 1,000th entry
    {
        double perc = count/(double)totalReadCount;
        std::lock_guard<std::mutex> guard(*printProgressLock);
        printProgress(perc);
    }

    return(result);
}

void MapAroundConstantBarcodesAsAnchorPolicy::map_pattern_between_linker(const std::string& seq, const int& oldEnd, 
                                                                         const int& start,
                                                                         BarcodePatternPtr barcodePatterns,
                                                                         std::vector<std::string>& barcodeList,
                                                                         int& barcodePosition, int& skippedBarcodes)
{
    //we have to map all barcodes that we missed
    int skipend=oldEnd;
    std::string skippedBarcodeString = seq.substr(skipend, start-skipend);
    int skipStringNewOffset = 0;
    for(int i = barcodePosition-skippedBarcodes; i < barcodePosition ; ++i)
    {
        BarcodeVector::iterator skippedPatternItr = barcodePatterns->begin() + i;
        std::string skippedBarcode = ""; //the actual real barcode that we find (mismatch corrected)
        int skipstart=0, skipscore = 0, skipdifferenceInBarcodeLength = 0;

        if(!(*skippedPatternItr)->match_pattern(skippedBarcodeString, 0, skipstart, skipend, skipscore, skippedBarcode, skipdifferenceInBarcodeLength, false, false, true))
        {
            barcodeList.push_back("");
        }
        else
        {
            skipStringNewOffset = skipend;
            skippedBarcodeString.erase(0,skipStringNewOffset);

            barcodeList.push_back(skippedBarcode);
        }
    }
}

bool MapAroundConstantBarcodesAsAnchorPolicy::split_line_into_barcode_patterns(const std::pair<fastqLine, fastqLine>& seq, 
                                        DemultiplexedLine& demultiplexedLine, const input& input, 
                                        BarcodePatternPtr barcodePatterns,
                                        fastqStats& stats)
{
    int offset = 0;
    int oldEnd = 0;

    std::vector<std::string> barcodeList;

    //firstly map each constant barcode
    int barcodePosition = 0;
    int skippedBarcodes = 0;
    for(BarcodeVector::iterator patternItr = barcodePatterns->begin(); 
        patternItr < barcodePatterns->end(); 
        ++patternItr)
    {

        //exclude non constant
        if(!(*patternItr)->is_constant())
        {
            ++barcodePosition;
            ++skippedBarcodes;
            continue;
        }
        
        //map next constant region
        int start=0, end=0, score = 0, differenceInBarcodeLength = 0;
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
            
        std::string subStringToSearchBarcodes = seq.first.line.substr(oldEnd, seq.first.line.length() - oldEnd);

        if(!(*patternItr)->match_pattern(subStringToSearchBarcodes, offset, start, end, score, barcode, differenceInBarcodeLength, false, false, true))
        {
            ++barcodePosition;
            ++skippedBarcodes;
            continue;
        }

        //if(start < oldEnd)
        //{
         //   ++stats.noMatches;
          //  return false;
        //} //we have in this case barcodes that are non sequential


        //write out all found barcodes
        //1.) write the skipped barcodes so far
        if(skippedBarcodes > 0)
        {
            map_pattern_between_linker(seq.first.line,oldEnd, start, barcodePatterns, demultiplexedLine.barcodeList, barcodePosition, skippedBarcodes);
            //std::string skippedBarcodeString = seq.first.substr(oldEnd, start-oldEnd);
            //barcodeList.push_back(skippedBarcodeString);
        }

        //2.) write the newly mapped constant barcode
        demultiplexedLine.barcodeList.push_back(barcode);

        //set paramters after constant mapping
        oldEnd += end;
        ++barcodePosition;
        skippedBarcodes = 0;
    }

    if(skippedBarcodes > 0)
    {
        int start = seq.first.line.length();
        map_pattern_between_linker(seq.first.line,oldEnd, start, barcodePatterns, demultiplexedLine.barcodeList, barcodePosition, skippedBarcodes);
        //std::string skippedBarcodeString = seq.first.substr(oldEnd, start-oldEnd);
        //barcodeList.push_back(skippedBarcodeString);
    }

    ++stats.perfectMatches;

    return true;
}

template class Mapping<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy>;
template class Mapping<MapEachBarcodeSequentiallyPolicyPairwise, ExtractLinesFromFastqFilePolicyPairedEnd>;
template class Mapping<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy>;
template class Mapping<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromTxtFilesPolicy>;
template class Mapping<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromFastqFilePolicy>;
