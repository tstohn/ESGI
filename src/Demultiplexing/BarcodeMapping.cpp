#include "BarcodeMapping.hpp"

bool check_if_seq_too_short(const size_t& offset, const std::string& seq)
{
    if(offset >= seq.length())
    {
        return true;
    }

    return false;
}

std::string_view trim(std::string_view sv)
{
    const char* begin = sv.data();
    const char* end   = sv.data() + sv.size();

    // Trim from the start
    while (begin < end && std::isspace(static_cast<unsigned char>(*begin))) {
        ++begin;
    }
    // Trim from the end
    while (end > begin && std::isspace(static_cast<unsigned char>(*(end - 1)))) {
        --end;
    }

    return std::string_view(begin, static_cast<size_t>(end - begin));
}

//parse a new file with barcodes and writes all barcodes into a vector
std::vector<std::string> parse_variable_barcode_file(const std::string& barcodeFile)
{
    std::vector<std::string> barcodes;
    try
    {
        std::ifstream barcodeFileStream;
        barcodeFileStream.open(barcodeFile);
        std::stringstream strStream;
        strStream << barcodeFileStream.rdbuf();
        std::string barcodeListString = strStream.str();
        //trim whitespaces (e.g., newline at end), we still check later for newlines at end of file, which is redundand for now    
        barcodeListString = trim(barcodeListString);
        if(barcodeListString.empty())
        {
            std::cerr << "Error: The abrcode file " << barcodeFile << " contains no barcodes. Please provide a list of comma seperated barcodes.\n";
        }

        size_t start = 0;
        size_t end = 0;
        while ((end = barcodeListString.find(',', start)) != std::string::npos) 
        {
            std::string_view seq(&barcodeListString[start], end - start);
            std::string_view trimmed = trim(seq);
            //validate string
            for (char const &c: trimmed) {
                if(!(c=='A' || c=='T' || c=='G' || c=='C' ||
                        c=='a' || c=='t' || c=='g' || c=='c'))
                        {   
                            if(c==' ' || c=='\t' || c=='\n')
                            {
                                std::cerr << "PARAMETER ERROR: Detected a whitespace in following barcode file: " << barcodeFile << ", pls. remove it to continue!\n";
                                std::cerr << "the whitespace occured within this barcode: " << trimmed << "\n";
                                std::cerr << "(Barcodes have to be comma-seperated with no whitespace, newline, tab in between. Also there should be no newline at the end of the file)\n";
                                if(c=='\n')
                                {
                                    std::cerr << "Barcodes should be comma seperated. Looks like you forgot to remove a newline between barcodes?\n";
                                }
                            }
                            else
                            {
                                std::cerr << "PARAMETER ERROR: a barcode sequence in " << barcodeFile << " is not a base (A,T,G,C)\n";
                            }
                            exit(1);
                        }
            }
            barcodes.emplace_back(trimmed);
            start = end + 1;
        }
        // last token
        std::string_view seq(&barcodeListString[start], barcodeListString.size() - start);
        std::string_view trimmed = trim(seq);
        // validate seq here
        for (char const &c: trimmed) 
        {
                if(!(c=='A' || c=='T' || c=='G' || c=='C' ||
                        c=='a' || c=='t' || c=='g' || c=='c'))
                        {   
                            if(c==' ' || c=='\t' || c=='\n')
                            {
                                std::cerr << "PARAMETER ERROR: Detected a whitespace in following barcode file: " << barcodeFile << ", pls. remove it to continue!\n";
                                std::cerr << "the whitespace occured within this barcode: " << trimmed << "\n";
                                std::cerr << "(Barcodes have to be comma-seperated with no whitespace, newline, tab in between. Also there should be no newline at the end of the file)\n";
                                if(c=='\n')
                                {
                                    std::cerr << "Barcodes should be comma seperated. Looks like you forgot to remove a newline between barcodes?\n";
                                }
                            }
                            else
                            {
                                std::cerr << "PARAMETER ERROR: a barcode sequence in " << barcodeFile << " is not a base (A,T,G,C)\n";
                            }
                            exit(1);
                        }
        }
        barcodes.emplace_back(trimmed);
        //close barcode file
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
std::pair<std::string, std::vector<std::string> >parse_pattern_line(std::string& line, int number)
{
    std::vector<std::string> patterns;
    std::string name;

    try {
        //if present parse the name of the line
        size_t colonPos = line.find(":");
        if (colonPos != std::string::npos) 
        {
            // Find the first '[' after the colon
            size_t bracketPos = line.find('[', colonPos);
            if (bracketPos != std::string::npos)
            {
                //make sure this patterns only exists once
                if (line.find(':', colonPos + 1) != std::string::npos) 
                {
                    std::cerr << "The substring <:> appears more than once. This is invalid as <:> should only be used to" <<
                    " name patterns like, e.g., <name:[GACGTA][pattern.txt]>... " << std::endl;
                    exit(EXIT_FAILURE);     
                }

                //parse name and trim it
                name = trim(line.substr(0, colonPos));
                line = trim(line.substr(bracketPos));  // Start from the first '['
            }
            else
            {
                std::cerr << "PARAMETER ERROR: Found ':' but no following pattern starting with '['\n";
                exit(1);
            }
        }
        else
        {
            name = "PATTERN_" + std::to_string(number);
        }
        name = '\'' + name + '\'';

        // parse one pattern line: e.g.: [ACGTTCAG]  [15X]   [file1.txt]
        std::string delimiter = "]";
        const char delimiter2 = '[';
        size_t pos = 0;
        std::string seq;
        //PARSE PATTERN
        while ((pos = line.find(delimiter)) != std::string::npos) {
            seq = trim(line.substr(0, pos));
            line = trim(line.erase(0, pos + 1));
            if(seq.empty() || seq.at(0) != delimiter2)
            {
                std::cerr << "PARAMETER ERROR: Wrong barcode pattern parameter, it looks like a \'[\' is missing\n";
                exit(1);
            }
            seq.erase(0, 1);  // remove the '['
            seq = trim(seq);  // trim the actual pattern content
            
            if (!seq.empty()) {  // only add non-empty patterns
                patterns.push_back(seq);
            }
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
std::vector<std::pair<std::string, std::vector<std::string>>> parse_pattern_file(const std::string& patternFile)
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
        std::string trimmedLine = trim(line);
        if(!line.empty())
        {
            patternList.emplace_back(parse_pattern_line(line, lineNumber));
            ++lineNumber;
        }
    }

    return(patternList);
}

//parse a new lnie in mismatch file
std::vector<std::vector<int>> parse_mismatch_file(const std::string& mismatchFile)
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
        //we trim every individual line in case whitespaces are present
        mismatchLine = trim(mismatchLine);
        if(mismatchLine.empty()){continue;}
        std::vector<int> mismatches;
        size_t pos = 0;
        while ((pos = mismatchLine.find(delimiter)) != std::string::npos) 
        {
            std::string seq = mismatchLine.substr(0, pos);
            mismatchLine.erase(0, pos + 1);
            if(seq.empty() || seq.find_first_not_of("0123456789") != std::string::npos)
            {
                std::cerr << "The mismatch file contains a non-numeric element [" << seq << "]\n";
                std::cerr << "Please provide comma seperated mismatches. One line for every barcode-pattern.\n";
                exit(EXIT_FAILURE);
            }
            mismatches.push_back(std::stoi(seq));
        }
        if(mismatchLine.empty() || mismatchLine.find_first_not_of("0123456789") != std::string::npos)
        {
            std::cerr << "The mismatch file contains a non-numeric element [" << mismatchLine << "]\n";
            std::cerr << "Please provide comma seperated mismatches. One line for every barcode-pattern.\n";
            exit(EXIT_FAILURE);
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
    std::unordered_map<std::string, std::vector<std::string>>& fileToBarcodesMap,
    const input& input)
{
    BarcodeVector barcodeVector;
    BarcodeVector independentReverseVector;
    bool isReversePattern = false;
    bool containsDNA = false;

    //iterate through all barcodes in the pattern line
    //parse the pattern and immediately create the Barcode
    for(size_t barcodeIdx = 0; barcodeIdx < barcodeList.size(); ++barcodeIdx)
    {

        std::string patternElement = barcodeList.at(barcodeIdx);
        int barcodeLength = 0;
        bool barcodeFound = false;

        std::cout << "  Creating barcodes for [" << patternElement << "]\n";
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
            if (input.independentReverseMapping && isReversePattern) 
            {
                independentReverseVector.push_back(barcodePtr);
            } else {
                barcodeVector.push_back(barcodePtr);
            }
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
            VariableBarcode barcode(fileToBarcodesMap.at(patternElement), patternElement, mismatchList.at(barcodeIdx), input.hamming);
            std::shared_ptr<VariableBarcode> barcodePtr(std::make_shared<VariableBarcode>(barcode));
            if (input.independentReverseMapping && isReversePattern) 
            {
                independentReverseVector.push_back(barcodePtr);
            } else {
                barcodeVector.push_back(barcodePtr);
            }
            barcodeFound = true;
        }

        //random element(e.g., UMI)
        //format: [15X], before we must check if the prefix-digit is not part of a file name (done above for variable barcode file)
        if(std::isdigit(patternElement[0]))
        {
            size_t pos = 0;
            barcodeLength = std::stoi(patternElement, &pos);
            std::string seq = patternElement.substr(pos);
            if(seq == "X" || seq == "x")
            {
                //For random sequences we MUST specific the length of this sequence
                WildcardBarcode barcode(mismatchList.at(barcodeIdx), patternElement, barcodeLength);
                std::shared_ptr<WildcardBarcode> barcodePtr(std::make_shared<WildcardBarcode>(barcode));
                if (input.independentReverseMapping && isReversePattern) 
                {
                    independentReverseVector.push_back(barcodePtr);
                } else 
                {
                    barcodeVector.push_back(barcodePtr);
                }
                barcodeFound = true;
            }
            else
            {
                //we already checked above if it is an existing barcode file
                std::cerr << "Encountered an error in barode pattern: " << patternElement << ".\n" <<
                "If a pattern has a digit as a prefix it must ether be part of a barcode-file name " <<
                "or it must be a random sequence like UMIs in the format <NUMBER><X> with a number followed by a single capital X." << 
                "For example [15X]\n";
                exit(1);
            }
        }

        //DNA element
        if(patternElement == "DNA")
        {
            DNABarcode barcode(mismatchList.at(barcodeIdx));
            std::shared_ptr<DNABarcode> barcodePtr(std::make_shared<DNABarcode>(barcode));
            if (input.independentReverseMapping && isReversePattern) 
            {
                independentReverseVector.push_back(barcodePtr);
            } else {
                barcodeVector.push_back(barcodePtr);
            }
            barcodeFound = true;
            containsDNA = true;
        }

        //stop pattern: [*] (only map up to here)
        if(patternElement == "*")
        {
            StopBarcode barcode(patternElement, mismatchList.at(barcodeIdx));
            std::shared_ptr<StopBarcode> barcodePtr(std::make_shared<StopBarcode>(barcode));
            if (input.independentReverseMapping && isReversePattern) 
            {
                independentReverseVector.push_back(barcodePtr);
            } else 
            {
                barcodeVector.push_back(barcodePtr);
            }
            barcodeFound = true;
        }

        //seperator pattern: [-] (ether stop here, or if DNA extract all till the read ends)
        if(patternElement == "-")
        {
            if (input.independentReverseMapping) 
            {
                isReversePattern = true;
            }
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
            "at this position (did you supply an existing path?), a random barcode like UMIs with <number><X> like <15X>, or a simple string such as <DNA> for a DNA/ RNA sequence," <<
            "or <*> for a stop position (to seperate FW, RV reads, or simply exclude a sequence in the middle).";
            exit(1);
        }
    }

    BarcodeVectorPtr barcodeVectorPtr = std::make_shared<BarcodeVector>(barcodeVector);
    BarcodePatternPtr patternPtr = std::make_shared<BarcodePattern>(BarcodePattern(containsDNA, patternName, barcodeVectorPtr));
    
    //set reverse pattern if present
    if (input.independentReverseMapping && !independentReverseVector.empty()) 
    {
        BarcodeVectorPtr reverseVectorPtr = std::make_shared<BarcodeVector>(independentReverseVector);
        patternPtr->independentReversePattern = reverseVectorPtr;
    }

    //if the independent flag is set but we never have a read-seperator throw an error
    if(input.independentReverseMapping && !isReversePattern)
    {
        std::cerr << "ERROR: To run in independent-mode (independent,seperate forward and reverse mapping, where we do not expect an overlap and therefore calcualte no reverse complement for the RV read) we need to give a [-] pattern:\n";
        std::cerr << "E.g.: PATTERN:[AAA][5X][-][TTT][5X]. This maps in the reverse read first TTT and then an UMI of 5 bases,\
        and it does not expect the reverse read to start with 5X.";
        exit(EXIT_FAILURE);
    }

    return(patternPtr);
}

template <typename MappingPolicy, typename FilePolicy>
bool Mapping<MappingPolicy, FilePolicy>::generate_barcode_patterns(const input& input)
{
    
    //a map of file names to the actual barcodes in this file. This map is updated while itereating through the pattern file
    //whenever an unknown file is encountered
    std::unordered_map<std::string, std::vector<std::string>> fileToBarcodesMap;

    //only parsing, no quality check
    std::vector<std::pair<std::string, std::vector<std::string>>> patternList = parse_pattern_file(input.barcodePatternsFile);
    std::vector<std::vector<int>> mismatchList = parse_mismatch_file(input.mismatchFile);

    // assert that the number of mismatches and patterns are the same
    if(patternList.size() != mismatchList.size())
    {
        std::cerr << "Number of barcode pattern lines and mismatch lines does not match. Please correct in parameters.\n";
        exit(1);
    }

    for(size_t i = 0; i < patternList.size(); i++)
    {
        if(patternList.at(i).second.size() != mismatchList.at(i).size())
        {
            std::cerr << "Number of barcode patterns and mismatches does not match for line " << std::to_string(i) << ". Please correct in parameters.\n";
            exit(EXIT_FAILURE);
        }
    }
    //make sure pattern names are unique
    std::unordered_map<std::string, int> countMap;
    for(size_t i = 0; i < patternList.size(); i++)
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
    for(size_t i = 0; i < patternList.size(); i++)
    {
        std::cout << "Creating pattern for " << patternList.at(i).first << "\n";
        barcodePatternList->emplace_back(create_barcodeVector_from_patternLine(patternList.at(i).second, mismatchList.at(i), patternList.at(i).first, fileToBarcodesMap, input));
    }

    return true;
}

//apply all BarcodePatterns to a single line to generate a vector strings (the Barcodemapping)
//realBarcodeMap contains the actual string in the seauence that gets a barcode assigned
bool MapEachBarcodeSequentiallyPolicy::split_line_into_barcode_patterns(const std::pair<fastqLine, fastqLine>& seq, 
                                        DemultiplexedLine& demultiplexedLine, const input& input, 
                                        BarcodePatternPtr barcodePatterns,
                                        int& mmScore,
                                        OneLineDemultiplexingStatsPtr stats)
{
    (void) input; //only needed for pairwise splitting of lines
    
    //fastq-read specific variables
    unsigned int positionInFastqLine = 0;
    int totalEdits = 0;
    int position = -1;

    for(BarcodeVector::iterator patternItr = barcodePatterns->begin(); 
        patternItr < barcodePatterns->end(); 
        ++patternItr)
    {

        //barcode-specific variables
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
        int del = 0; 
        int ins = 0;
        int subst = 0;
        int targetEnd = 0;
        ++position;

        //if we have a wildcard, just extract the necessary barcodes
        if((*patternItr)->is_wildcard())
        {
            (*patternItr)->align(barcode, seq.first.line, positionInFastqLine, targetEnd, del, ins, subst);
            positionInFastqLine += (targetEnd);
            demultiplexedLine.barcodeList.push_back(barcode);

            if(stats != nullptr)
            {
                stats->insertions.push_back(0);
                stats->deletions.push_back(0);
                stats->substitutions.push_back(0);
            }

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
            demultiplexedLine.dna = seq.first.line.substr(positionInFastqLine, seq.first.line.length()); //extract DNA fragment
            demultiplexedLine.dnaQuality = seq.first.quality.substr(positionInFastqLine, seq.first.line.length()); //extract DNA fragment
            continue;
        }
        else if((*patternItr)->is_read_end())
        {
            //stop here: we do not continue mapping when the read ends
            break;
        }

        // could happen in the case of deletions in the UMI sequence...
        bool seqToShort = check_if_seq_too_short(positionInFastqLine, seq.first.line);
        if(seqToShort)
        {
            //std::cout << "too short: " << positionInFastqLine << " " << seq.first.line.length() << "\n";
            //++stats->noMatches;
            return false;
        }

        //std::cout << " trying " << seq.first.line << "\n";
        if(!(*patternItr)->align(barcode, seq.first.line, positionInFastqLine, targetEnd, del, ins, subst))
        {            
            //save until where we mapped
            if(stats != nullptr)
            {
                stats->failedLinesMappingFw.first = barcodePatterns->patternName;
                stats->failedLinesMappingFw.second = position;
            }
            return false;
        }

        //std::cout << "ALIGN: found " << barcode << " within " << seq.first.line.substr(positionInFastqLine,targetEnd) << " with D: " << del << ", I: " << ins << " and S: " << subst << "\n";

        totalEdits = totalEdits + del + ins + subst;
        positionInFastqLine += targetEnd; //targetEnd is zero indexed alst position in target-sequence that maps to pattern
        //positionInFastqLine is the first position to INCLUDE in next alignment

        assert(barcode != "");
        if(stats != nullptr)
        {
            stats->insertions.push_back(ins);
            stats->deletions.push_back(del);
            stats->substitutions.push_back(subst);
        }
        
        //add this match to the BarcodeMapping
        //barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        demultiplexedLine.barcodeList.push_back(barcode);
    }

    if(stats != nullptr)
    {

        if(totalEdits == 0)
        {
            ++stats->perfectMatches;
        }
        else
        {
            ++stats->moderateMatches;
        }
    }
    mmScore = totalEdits;

    return true;
}

bool MapEachBarcodeSequentiallyPolicyPairwise::map_forward(const fastqLine& seq, 
                                                           BarcodePatternPtr barcodePatterns,
                                                           OneLineDemultiplexingStatsPtr stats,
                                                           DemultiplexedLine& demultiplexedLine,
                                                           unsigned int& barcodePosition,
                                                           int& totalEdits,
                                                           PatternType type = PatternType::Forward)
{

    //fastq-read specific variables
    unsigned int positionInFastqLine = 0;
    int position = -1;

    for(BarcodeVector::iterator patternItr = barcodePatterns->begin(type); 
        patternItr < barcodePatterns->end(type); 
        ++patternItr)
    {

        //barcode-specific variables
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
        int del = 0; 
        int ins = 0;
        int subst = 0;
        int targetEnd = 0;
        ++position;

        //if we have mapped to the end of the sequence (but barcodes of the pattern are still missing)
        if( (seq.line.size() <= positionInFastqLine) && !(*patternItr)->is_stop() && !(*patternItr)->is_read_end()){return false;}

        //if we have a wildcard skip this matching, we match again the next sequence
        if((*patternItr)->is_wildcard())
        {
            (*patternItr)->align(barcode, seq.line, positionInFastqLine, targetEnd, del, ins, subst);
            positionInFastqLine += (targetEnd);
            ++barcodePosition; //increase the count of found positions
            demultiplexedLine.barcodeList.push_back(barcode);

            if(stats != nullptr)
            {       
                stats->insertions.push_back(0);
                stats->deletions.push_back(0);
                stats->substitutions.push_back(0);
            }
            continue;
        }
        else if((*patternItr)->is_stop())
        {

            //the stop positions is also a found position (count only for fw - otherwise sequences could wrongfully overlap)
            ++barcodePosition; //increase the count of found positions
            //stop here: we do not continue mapping after stop barcode [*]
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
            demultiplexedLine.dna = seq.line.substr(positionInFastqLine, seq.line.length()); //extract DNA fragment
            demultiplexedLine.dnaQuality = seq.quality.substr(positionInFastqLine, seq.line.length()); //extract DNA fragment
            continue;
        }
        else if((*patternItr)->is_read_end())
        {
            ++barcodePosition; //read-end is a valid found barcode (count only for forward - otherwise sequences could wrongfully overlap)
            //stop here: we do not continue mapping when the read ends
            return true;
        }

        // in case the sequence ends in the middle of the next barcode
        bool seqToShort = check_if_seq_too_short(positionInFastqLine, seq.line);
        if(seqToShort)
        {
            return false;
        }

        //IF NON OF THE ABOVE - TRY TO MAP PATTERN

        //if we did not match a pattern
        if(!(*patternItr)->align(barcode,seq.line, positionInFastqLine, targetEnd, del, ins, subst))
        {
            //save until where we mapped
            if(stats != nullptr)
            {
                stats->failedLinesMappingFw.first = barcodePatterns->patternName;
                stats->failedLinesMappingFw.second = position;
            }
            return false;
        }
        
        //std::cout << "ALIGN FW: found " << barcode << " at start " << positionInFastqLine << " with new end " << targetEnd << "\n";

        totalEdits = totalEdits + del + ins + subst;
        positionInFastqLine += targetEnd; //targetEnd is zero indexed alst position in target-sequence that maps to pattern
        //positionInFastqLine is the first position to INCLUDE in next alignment

        assert(barcode != "");
        if(stats != nullptr)
        {
            stats->insertions.push_back(ins);
            stats->deletions.push_back(del);
            stats->substitutions.push_back(subst);
        }

        //add this match to the BarcodeMapping
        //barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        demultiplexedLine.barcodeList.push_back(barcode);
        ++barcodePosition; //increase the count of found positions
    }

    return true;
}

bool MapEachBarcodeSequentiallyPolicyPairwise::map_reverse(const fastqLine& seq, 
                                                           BarcodePatternPtr barcodePatterns,
                                                           OneLineDemultiplexingStatsPtr stats,
                                                           DemultiplexedLine& demultiplexedLine,
                                                           unsigned int& barcodePosition,
                                                           int& totalEdits)
{
    //fastq-read specific variables
    unsigned int positionInFastqLine = 0;
    int position = -1;

    //iterate reverse through patterns
    for(BarcodeVector::reverse_iterator patternItr = barcodePatterns->rbegin(); 
        patternItr < barcodePatterns->rend(); 
        ++patternItr)
    {

        //barcode-specific variables
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
        int del = 0; 
        int ins = 0;
        int subst = 0;
        int targetEnd = 0;
        ++position;

        //if we have mapped to the end of the sequence (but barcodes of the pattern are still missing)
        if( (seq.line.size() <= positionInFastqLine) && !(*patternItr)->is_stop() && !(*patternItr)->is_read_end()){return false;}

        //if we have a wildcard skip this matching, we match again the next sequence
        if((*patternItr)->is_wildcard())
        {
            (*patternItr)->align(barcode, seq.line, positionInFastqLine, targetEnd, del, ins, subst);
            positionInFastqLine += (targetEnd);
            ++barcodePosition; //increase the count of found positions
            demultiplexedLine.barcodeList.push_back(barcode);

            if(stats != nullptr)
            {
                stats->insertions.push_back(0);
                stats->deletions.push_back(0);
                stats->substitutions.push_back(0);
            }
            continue;
        }
        else if((*patternItr)->is_stop())
        {
            //stop here: we do not continue mapping after stop barcode [*]
            // DO NOT increase the count of found positions (++barcodePosition;),
            //it is ONLY counted in the forward read
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
            demultiplexedLine.dna = seq.line.substr(positionInFastqLine, seq.line.length()); //extract DNA fragment
            demultiplexedLine.dnaQuality = seq.quality.substr(positionInFastqLine, seq.line.length()); //extract DNA fragment
            continue;
        }
        else if((*patternItr)->is_read_end())
        {
            // DO NOT increase the count of found positions (++barcodePosition;), 
            //it is ONLY counted in the forward read
            //stop here: we do not continue mapping when the read ends
            return true;
        }

        // in case the sequence ends in the middle of the next barcode
        bool seqToShort = check_if_seq_too_short(positionInFastqLine, seq.line);
        if(seqToShort)
        {
            return false;
        }

        //map each pattern with reverse complement
        if(!(*patternItr)->align(barcode,seq.line, positionInFastqLine, targetEnd, del, ins, subst, true))
        {
            //save until where we mapped
            if(stats != nullptr)
            {
                stats->failedLinesMappingRv.first = barcodePatterns->patternName;
                stats->failedLinesMappingRv.second = position;
            }
            return false;
        }
        //std::cout << "ALIGN RV: found " << barcode << " at start " << positionInFastqLine << " with new end " << targetEnd << "\n";

        totalEdits = totalEdits + del + ins + subst;
        positionInFastqLine += targetEnd; //targetEnd is zero indexed alst position in target-sequence that maps to pattern
        //positionInFastqLine is the first position to INCLUDE in next alignment

        assert(barcode != "");
        if(stats != nullptr)
        {
            stats->insertions.push_back(ins);
            stats->deletions.push_back(del);
            stats->substitutions.push_back(subst);
        }
      
        //add this match to the BarcodeMapping
        //barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        demultiplexedLine.barcodeList.push_back(barcode);
        ++barcodePosition; //increase the count of found positions
    }

    return true;
}

//demultiplexedLineFw is actually already the final barcodeList (stored in Demultiplexline), to not copy vectors of strings around but add in place
//if reverse read was DNA and contains that information we need to copy it into demultiplexedLineFw
bool MapEachBarcodeSequentiallyPolicyPairwise::combine_mapping(const BarcodePatternPtr& barcodePatterns,
                                                               DemultiplexedLine& demultiplexedLineFw,
                                                               const unsigned int& barcodePositionFw,
                                                               const DemultiplexedLine& demultiplexedLineRv,
                                                               const unsigned int& barcodePositionRv,
                                                               OneLineDemultiplexingStatsPtr stats,
                                                               OneLineDemultiplexingStatsPtr statsRv,
                                                               int& score_sum)
{
    //check that we span the whole sequence (except constant regions)
    unsigned int patternNum = barcodePatterns->size();

    //if the reverse read contains DNA co-y the information into the Fw demultiplexed line (the fw line is used as final class)
    if(demultiplexedLineRv.dna != "")
    {
        demultiplexedLineFw.dna = demultiplexedLineRv.dna; //extract DNA fragment
        demultiplexedLineFw.dnaQuality = demultiplexedLineRv.dnaQuality; //extract DNA fragment
    }

    //if positions are next to each other just return
    //in case of a stop pattern [*], we added +1 to the barcodePositionFw, so that barcodePositionFw+barcodePositionRv should be euqual to patternNum
    if(patternNum == (barcodePositionFw + barcodePositionRv))
    {
        for (int i = static_cast<int>(demultiplexedLineRv.barcodeList.size()) - 1; i >= 0; --i)
        {
            demultiplexedLineFw.barcodeList.push_back({demultiplexedLineRv.barcodeList.at(i)});
            if(stats != nullptr)
            {
                stats->insertions.push_back(statsRv->insertions.at(i));
                stats->deletions.push_back(statsRv->deletions.at(i));
                stats->substitutions.push_back(statsRv->substitutions.at(i));

                //reset values for failed barcode mapping positions if we mapped perfectly
                //it could be that in, e.g. the forward read we did not map a last barcode but it was mapped in
                // the reverse read...
                stats->failedLinesMappingFw.first = "";
                stats->failedLinesMappingFw.second = 0;
                stats->failedLinesMappingRv.first = "";
                stats->failedLinesMappingRv.second = 0;
            }
        }
    }

    //if they r overlapping, check all overlapping ones are the same
    if(patternNum < (barcodePositionFw + barcodePositionRv))
    {
        //assert all overlapping barcodes r the same
        size_t start = (patternNum - barcodePositionRv); // this is the idx of the first shared barcode
        size_t end = barcodePositionFw - 1; //this is the idnex of the last shared barcode
        size_t j = demultiplexedLineRv.barcodeList.size() - 1;
        for(size_t i = start; i <= end; ++i)
        {
            if(demultiplexedLineFw.barcodeList.at(i) != demultiplexedLineRv.barcodeList.at(j))
            {
                //++stats->noMatches;
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
            if(stats != nullptr)
            {
                stats->insertions.push_back(statsRv->insertions.at(i));
                stats->deletions.push_back(statsRv->deletions.at(i));
                stats->substitutions.push_back(statsRv->substitutions.at(i));
            }
        }
    }

    //if there r barcodePatterns missing in the middle, check they r only constant and add any
    //non-found-constant barcode
    if(patternNum > (barcodePositionFw + barcodePositionRv))
    {
        //check if missing barcodes are only constant
        size_t start = barcodePositionFw; // index of first missing barcode
        size_t end = (patternNum - barcodePositionRv) - 1; //last idnex that can be missing

        for(size_t i = start; i <= end; ++i)
        {
            if(!(barcodePatterns->barcodePattern->at(i)->is_constant()))
            {
                //++stats->noMatches;
                return false;
            }
            demultiplexedLineFw.barcodeList.push_back(*(barcodePatterns->barcodePattern->at(i)->get_patterns().at(0)));

            //in the statistics save this as zero MM, since we can not accurately count the number of MM
            if(stats != nullptr)
            {
                stats->insertions.push_back(0);
                stats->deletions.push_back(0);
                stats->substitutions.push_back(0);
            }
        }

        //add reverse patterns        
        for (int i = static_cast<int>(demultiplexedLineRv.barcodeList.size()) - 1; i >= 0; --i)
        {
            demultiplexedLineFw.barcodeList.push_back({demultiplexedLineRv.barcodeList.at(i)});
            if(stats != nullptr)
            {
                stats->insertions.push_back(statsRv->insertions.at(i));
                stats->deletions.push_back(statsRv->deletions.at(i));
                stats->substitutions.push_back(statsRv->substitutions.at(i));
            }
        }

        //reset values for failed barcode mapping positions, we assume we correctly mapped this one
        if(stats != nullptr)
        {
            stats->failedLinesMappingFw.first = "";
            stats->failedLinesMappingFw.second = 0;
            stats->failedLinesMappingRv.first = "";
            stats->failedLinesMappingRv.second = 0;
        }
    }

    if(stats != nullptr)
    {
        if(score_sum == 0)
        {
            ++stats->perfectMatches;
        }
        else
        {
            ++stats->moderateMatches;
        }
    }

    return true;
}

bool MapEachBarcodeSequentiallyPolicyPairwise::split_line_into_barcode_patterns(const std::pair<fastqLine, fastqLine>& seq, 
                                        DemultiplexedLine& demultiplexedLine,
                                        const input& input,
                                        BarcodePatternPtr barcodePatterns,
                                        int& mmScore,
                                        OneLineDemultiplexingStatsPtr stats)
{

    bool pairwiseMappingSuccess = false;

    //if both foward and reverse read are mapped from 5' - 3' and reverse read is no reverse complement
    if(input.independentReverseMapping)
    {

        unsigned int barcodePositionFw = 0;
        int tmpMMScore = 0;
        //for forward read we immediately add barcodes to demultiplexedLine.barcodeList, which is then extended in combine pattern, IF we find all patterns
        bool forwardSuccess = map_forward(seq.first, barcodePatterns, stats, demultiplexedLine, barcodePositionFw, tmpMMScore);

        DemultiplexedLine demultiplexedLineRv;
        unsigned int barcodePositionRv = 0;

        //statistics result for the reverse line
        //the forward line is saved in the final result stat object, this object is then later
        //updated with the reverse read information
        OneLineDemultiplexingStatsPtr statsRvPtr;
        //if we save stats initialize the temporary one here
        if(stats != nullptr){statsRvPtr = std::make_shared<OneLineDemultiplexingStats>();}
        else{statsRvPtr = nullptr;}

        bool reverseSuccess = map_forward(seq.second, barcodePatterns, statsRvPtr, demultiplexedLineRv,barcodePositionRv, tmpMMScore, PatternType::Reverse);

        if(forwardSuccess && reverseSuccess)
        {
            mmScore = tmpMMScore;
            pairwiseMappingSuccess = true;
            for (int i =  0; i < static_cast<int>(demultiplexedLineRv.barcodeList.size()); ++i)
            {
                demultiplexedLine.barcodeList.push_back({demultiplexedLineRv.barcodeList.at(i)});
                if(stats != nullptr)
                {
                    stats->insertions.push_back(statsRvPtr->insertions.at(i));
                    stats->deletions.push_back(statsRvPtr->deletions.at(i));
                    stats->substitutions.push_back(statsRvPtr->substitutions.at(i));

                    //reset values for failed barcode mapping positions if we mapped perfectly
                    //it could be that in, e.g. the forward read we did not map a last barcode but it was mapped in
                    // the reverse read...
                    stats->failedLinesMappingFw.first = "";
                    stats->failedLinesMappingFw.second = 0;
                    stats->failedLinesMappingRv.first = "";
                    stats->failedLinesMappingRv.second = 0;
                }
            }

            if(stats != nullptr)
            {
                if(mmScore == 0)
                {
                    ++stats->perfectMatches;
                }
                else
                {
                    ++stats->moderateMatches;
                }
            }
        }
    }
    else 
    {
        //map forward reads barcodeList contains the stored barcodes, 
        //barcodePosition is the psoiton of the last mapped barcode
        unsigned int barcodePositionFw = 0;
        int tmpMMScore = 0;
        //for forward read we immediately add barcodes to demultiplexedLine.barcodeList, which is then extended in combine pattern, IF we find all patterns
        map_forward(seq.first, barcodePatterns, stats, demultiplexedLine, barcodePositionFw, tmpMMScore);
        DemultiplexedLine demultiplexedLineRv;
        unsigned int barcodePositionRv = 0;

    //  std::cout << "FOUND FOWARD: ";
    //  for(auto el : demultiplexedLine.barcodeList)
    //  {
    //      std::cout << el << "\t";
    //  }
    //  std::cout << "\n";

        //statistics result for the reverse line
        //the forward line is saved in the final result stat object, this object is then later
        //updated with the reverse read information
        OneLineDemultiplexingStatsPtr statsRvPtr;
        //if we save stats initialize the temporary one here
        if(stats != nullptr){statsRvPtr = std::make_shared<OneLineDemultiplexingStats>();}
        else{statsRvPtr = nullptr;}
        
        map_reverse(seq.second, barcodePatterns, statsRvPtr, demultiplexedLineRv,barcodePositionRv, tmpMMScore);

    // std::cout << "FOUND REVERSE: ";
    // for(auto el : demultiplexedLineRv.barcodeList)
    // {
    //     std::cout << el << "\t";
    // }
    // std::cout << "\n";

        //todo: better check if fw and rv mapped: at the moment they returna  fail if e.g. a half-constant barcode does not map
        pairwiseMappingSuccess = combine_mapping(barcodePatterns, demultiplexedLine, barcodePositionFw, demultiplexedLineRv, barcodePositionRv, stats, statsRvPtr, tmpMMScore);
        if(pairwiseMappingSuccess){mmScore = tmpMMScore;}
    }

    return (pairwiseMappingSuccess);
}

template <typename MappingPolicy, typename FilePolicy>
bool Mapping<MappingPolicy, FilePolicy>::demultiplex_read(const std::pair<fastqLine, fastqLine>& seq, 
                                                          DemultiplexedLine& demultiplexedLine,
                                                          BarcodePatternPtr pattern,
                                                          const input& input, 
                                                          std::atomic<unsigned long long>& count, const unsigned long long& totalReadCount,
                                                          int& mmScore,
                                                          OneLineDemultiplexingStatsPtr stats)
{
    //split line into patterns (barcodeMap, barcodePatters, stats are passed as reference or ptr)
    //and can be read by each thread, "addValue" method for barcodeMap is thread safe also for concurrent writing
    bool result;

    //demultipelxed barcodes are stored in barcodeMap
    result = this->split_line_into_barcode_patterns(seq, demultiplexedLine, input, pattern, mmScore, stats);

    //update status bar
    unsigned long long step = (totalReadCount + 50) / 100;
    if (step > 0 && count % step == 0 && totalReadCount != ULLONG_MAX && count < totalReadCount) //update at every 1,000th entry
    {
        double perc = count/(double)totalReadCount;
        std::lock_guard<std::mutex> guard(*printProgressLock);
        printProgress(perc);
    }


    return(result);
}


//THIS POLICY IS FOR NOW DISABLED, IT NEEEDS MAYOR CHANGES AFTER USING NEW ALIGNER AND IS NO LONGER SUPPORTED
/*
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
    for(size_t i = barcodePosition-skippedBarcodes; i < barcodePosition ; ++i)
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
                                        int& mmScore,
                                        OneLineDemultiplexingStatsPtr stats)
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

        //set parameters after constant mapping
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

    if(stats != nullptr)
    {
        ++stats->perfectMatches;
    }

    return true;
}
*/

template class Mapping<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy>;
template class Mapping<MapEachBarcodeSequentiallyPolicyPairwise, ExtractLinesFromFastqFilePolicyPairedEnd>;
template class Mapping<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy>;

//ANCHOR POLICIES ARE NOT LONGER SUPPORTED, LAST WORKING COMMIT IS BEFORE MERGE OF MULTIPATTERN BRANCH
//template class Mapping<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromTxtFilesPolicy>;
//template class Mapping<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromFastqFilePolicy>;
