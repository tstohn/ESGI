
#include "BarcodeMapping.hpp"

bool check_if_seq_too_short(const int& offset, const std::string& seq)
{
    if(offset >= seq.length())
    {
        return true;
    }

    return false;
}

//initialize the dictionary of mismatches in each specific barcode
// the vector in the dict has length "mismatches + 2" for each barcode
// one entry for zero mismatches, eveery number from, 1 to mismatches and 
//one for more mismatches than the max allowed number of mismathces
template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::initializeStats()
{
    for(BarcodePatternVector::iterator patternItr = barcodePatterns->begin(); 
        patternItr < barcodePatterns->end(); 
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

template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::parse_barcode_data(const input& input, std::vector<std::pair<std::string, char> >& patterns, 
                                                            std::vector<int>& mismatches, std::vector<std::vector<std::string> >& varyingBarcodes)
{
    try{
        // parse the pattern, mismatches, and barcode file (perform quality check as well)
        int numberOfNonConstantBarcodes = 0;
        std::string pattern = input.patternLine;
        std::string delimiter = "]";
        const char delimiter2 = '[';
        size_t pos = 0;
        std::string seq;
        //PARSE PATTERN
        std::pair<std::string, bool> seqPair;
        while ((pos = pattern.find(delimiter)) != std::string::npos) {
            seq = pattern.substr(0, pos);
            pattern.erase(0, pos + 1);
            if(seq.at(0) != delimiter2)
            {
                std::cerr << "PARAMETER ERROR: Wrong barcode patter parameter, it looks like a \'[\' is missing\n";
                exit(1);
            }
            seq.erase(0, 1);
            //just check if we have a barcode pattern of only'N', bcs the nunber of those patterns must match the number of lines in barcodefile
            bool nonConstantSeq = true;
            char patternType = 'c';
            for (char const &c: seq) {
                if(!(c=='A' || c=='T' || c=='G' || c=='C' ||
                    c=='a' || c=='t' || c=='g' || c=='c' ||
                    c=='N' || c=='X'))
                {
                    std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C,N)\n";
                    if(c==' ' || c=='\t' || c=='\n')
                    {
                        std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
                    }
                    exit(1);
                }
                if(c!='N'){nonConstantSeq=false;}

                //determine pattern type and throw error
                if( (patternType!='c') && (c!='N') && (c!='X') )
                {
                    std::cerr << "PARAMETER ERROR: a barcode sequence has bases as well as wildcard 'X' or variable 'N' bases, the combination is not allowed!!!\n";
                    exit(1);
                }
                if(c=='N'){patternType='v';}
                else if(c=='X'){patternType='w';}
                else if(c=='D'){patternType='d';}
            }
            if(nonConstantSeq){++numberOfNonConstantBarcodes;}
            
            patterns.push_back(std::make_pair(seq, patternType));
        }
        //PARSE mismatches
        pattern = input.mismatchLine;
        delimiter = ",";
        pos = 0;
        while ((pos = pattern.find(delimiter)) != std::string::npos) {
            seq = pattern.substr(0, pos);
            pattern.erase(0, pos + 1);
            mismatches.push_back(stoi(seq));
        }
        mismatches.push_back(stoi(pattern));
        if(patterns.size() != mismatches.size())
        {
            std::cerr << "PARAMETER ERROR: Number of barcode patterns and mismatches is not equal\n";
            exit(1);
        }
        //PARSE barcode file
        std::ifstream barcodeFile(input.barcodeFile);
        for(std::string line; std::getline(barcodeFile, line);)
        {
            delimiter = ",";
            pos = 0;
            std::vector<std::string> seqVector;
            while ((pos = line.find(delimiter)) != std::string::npos) {
                seq = line.substr(0, pos);
                line.erase(0, pos + 1);
                for (char const &c: seq) {
                    if(!(c=='A' | c=='T' | c=='G' |c=='C' |
                         c=='a' | c=='t' | c=='g' | c=='c'))
                         {
                            std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C)\n";
                            if(c==' ' | c=='\t' | c=='\n')
                            {
                                std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
                            }
                            exit(1);
                         }
                }
                seqVector.push_back(seq);
            }
            seq = line;
            for (char const &c: seq) {
                if(!(c=='A' || c=='T' || c=='G' || c=='C' ||
                        c=='a' || c=='t' || c=='g' || c=='c'))
                        {
                        std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C)\n";
                        if(c==' ' || c=='\t' || c=='\n')
                        {
                            std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence of barcode file; remove it to continue!\n";
                        }
                        exit(1);
                        }
            }
            seqVector.push_back(seq);
            varyingBarcodes.push_back(seqVector);
            seqVector.clear();
        }
        barcodeFile.close();
        if(numberOfNonConstantBarcodes != varyingBarcodes.size())
        {
            std::cerr << "PARAMETER ERROR: Number of barcode patterns for non-constant sequences [N*] and lines in barcode file are not equal\n";
            std::cerr << std::to_string(numberOfNonConstantBarcodes) << " and " << std::to_string(varyingBarcodes.size()) << " respectively\n";
            exit(1);
        }
        //PARSE guide file (only if we map AB and guide reads simultaneously)
        if(input.guideFile!="")
        {
            std::ifstream guideFile(input.guideFile);
            int lineCount = 0;
            std::vector<std::string> seqVector;
            for(std::string line; std::getline(guideFile, line);)
            {
                delimiter = ",";
                pos = 0;
                while ((pos = line.find(delimiter)) != std::string::npos) {
                    //if we find delimiters in more than one line, the file is not correct
                    if(lineCount > 0)
                    {
                        std::cerr << "PARAMETER ERROR: The input guide file contains too many lines. It should contain a single line with comma seperated guide sequences.\n";
                        exit(1);
                    }
                    seq = line.substr(0, pos);
                    line.erase(0, pos + 1);
                    for (char const &c: seq) {
                        if(!(c=='A' | c=='T' | c=='G' |c=='C' |
                            c=='a' | c=='t' | c=='g' | c=='c'))
                            {
                                std::cerr << "PARAMETER ERROR: a barcode sequence in guide file is not a base (A,T,G,C)\n";
                                if(c==' ' | c=='\t' | c=='\n')
                                {
                                    std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence of guide file; remove it to continue!\n";
                                }
                                exit(1);
                            }
                    }
                    seqVector.push_back(seq);
                }
                seq = line;
                for (char const &c: seq) {
                    if(!(c=='A' || c=='T' || c=='G' || c=='C' ||
                            c=='a' || c=='t' || c=='g' || c=='c'))
                            {
                            std::cerr << "PARAMETER ERROR: a barcode sequence in guide file is not a base (A,T,G,C)\n";
                            if(c==' ' || c=='\t' || c=='\n')
                            {
                                std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence of guide file; remove it to continue!\n";
                            }
                            exit(1);
                            }
                }
                seqVector.push_back(seq);
                ++lineCount;
            }
            guideList = std::make_shared<std::vector<std::string>>(seqVector);
            guideFile.close();
            //also set the bool storing weather guide reads contain a UMI or not
            guideUMI = input.guideUMI;
        }
    }
    catch(std::exception& e)
    {
        std::cerr << "PARAMETER ERROR: Check the input parameters again (barcode pattern list, mismatch list, file with barcode lists)\n";
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

template <typename MappingPolicy, typename FilePolicy>
std::vector<std::pair<std::string, char> > Mapping<MappingPolicy, FilePolicy>::generate_barcode_patterns(const input& input)
{
    //temporary structure storing the order, length and type of each pattern
    std::vector<std::pair<std::string, char> > patterns; // vector of all string patterns, 
                                                        //second entry is c=constant, v=varying, w=wildcard

    std::vector<int> mismatches; // vector of all string patterns
    std::vector<std::vector<std::string> > varyingBarcodes; // a vector storing for all non-constant barcode patterns in the order of occurence
                                                            // in the barcode pattern string the possible barcode sequences
    //fill the three vectors and handle as many errors as possible
    parse_barcode_data(input, patterns, mismatches, varyingBarcodes);

    //iterate over patterns and fill BarcodePatternVector instance
    BarcodePatternVector barcodeVector;
    BarcodePatternVector guideBarcodeVector;
    int variableBarcodeIdx = 0; // index for variable barcode sequences is shorter than the vector of patterns in total
    for(int i=0; i < patterns.size(); ++i)
    {
        if(patterns.at(i).second=='v')
        {
            VariableBarcode barcode(varyingBarcodes.at(variableBarcodeIdx), mismatches.at(i));
            std::shared_ptr<VariableBarcode> barcodePtr(std::make_shared<VariableBarcode>(barcode));
            barcodeVector.push_back(barcodePtr);
            
            //if we are at the AB/GUIDE barcode position, we need to hear insert the guide sequences
            if(variableBarcodeIdx == input.guidePos)
            {
                VariableBarcode barcode(*guideList, mismatches.at(i));
                std::shared_ptr<VariableBarcode> barcodePtr(std::make_shared<VariableBarcode>(barcode));
                guideBarcodeVector.push_back(barcodePtr);
            }
            else
            {
                guideBarcodeVector.push_back(barcodePtr);
            }
            ++variableBarcodeIdx;
        }
        else if(patterns.at(i).second=='c')
        {
            ConstantBarcode barcode(patterns.at(i).first, mismatches.at(i));
            std::shared_ptr<ConstantBarcode> barcodePtr(std::make_shared<ConstantBarcode>(barcode));
            barcodeVector.push_back(barcodePtr);
            guideBarcodeVector.push_back(barcodePtr);
        }
        else if(patterns.at(i).second=='w')
        {
            WildcardBarcode barcode(patterns.at(i).first, mismatches.at(i));
            std::shared_ptr<WildcardBarcode> barcodePtr(std::make_shared<WildcardBarcode>(barcode));
            barcodeVector.push_back(barcodePtr);
            //only add the wildcard to guidebarcodePattern if guide has a UMI 
            //(assuming the wildcard is only the UMI in this scenario)
            if(guideUMI)
            {
                guideBarcodeVector.push_back(barcodePtr);
            }
        }
        else if(patterns.at(i).second=='d')
        {
            //ERROR NOT IMPLEMENTED YET
        }

    }
    BarcodePatternVectorPtr barcodePatternVector = std::make_shared<BarcodePatternVector>(barcodeVector);

    //check that number of mismatches does not exceed number of bases in patter
    for(int i = 0; i < barcodePatternVector->size(); i++)
    {
        if(barcodePatternVector->at(i)->get_patterns().at(0).length() <= barcodePatternVector->at(i)->mismatches)
        {
            std::cerr << "Number of mismatches should not exceed the number of bases in pattern. Please correct in parameters.\n";
            exit(1);
        }
    }

    //set the vector of barcode patterns
    barcodePatterns = barcodePatternVector;
    guideBarcodePatterns = std::make_shared<BarcodePatternVector>(guideBarcodeVector);
    return patterns;
}

//apply all BarcodePatterns to a single line to generate a vector strings (the Barcodemapping)
//realBarcodeMap contains the actual string in the seauence that gets a barcode assigned
bool MapEachBarcodeSequentiallyPolicy::split_line_into_barcode_patterns(std::pair<const std::string&, const std::string&> seq, const input& input, 
                                        DemultiplexedReads& barcodeMap, 
                                        BarcodePatternVectorPtr barcodePatterns,
                                        fastqStats& stats)
{
    std::vector<std::string> barcodeList;

    //iterate over BarcodeMappingVector
    int offset = 0;
    int score_sum = 0;

    int old_offset = offset;
    unsigned int wildCardToFill = 0;
    int wildCardLength = 0, differenceInBarcodeLength = 0;
    for(BarcodePatternVector::iterator patternItr = barcodePatterns->begin(); 
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
            wildCardLength = (*patternItr)->get_patterns().at(0).length();
            offset += wildCardLength;
            wildCardToFill += 1;
            continue;
        }
        //for every barcodeMapping element find a match
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
        int start=0, end=0, score = 0;
        bool startCorrection = false;
        if(wildCardToFill){startCorrection = true;} // sart correction checks if we have to move our mapping window to the 5' direction
        // could happen in the case of deletions in the UMI sequence...
        bool seqToShort = check_if_seq_too_short(offset, seq.first);
        if(seqToShort)
        {
            return false;
        }
        if(!(*patternItr)->match_pattern(seq.first, offset, start, end, score, barcode, differenceInBarcodeLength, startCorrection, false))
        {
            ++stats.noMatches;
            return false;
        }
        
        std::string sequence = seq.first.substr(offset + start, end-start);
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
            std::string oldWildcardMappedBarcode = seq.first.substr(old_offset, (offset+start-end) - old_offset);
            //barcodeMap.emplace_back(uniqueChars.  (oldWildcardMappedBarcode));
            unsigned int lastWildcardEnd = 0; //in case we have several wildcards and need to know old offset
            // to subset new string
            while(wildCardToFill != 0)
            {
                BarcodePatternVector::iterator wildCardIt = patternItr;
                int pos = -1 * wildCardToFill;
                std::advance(wildCardIt, pos);
                int currentWildcardLength = (*wildCardIt)->get_patterns().at(0).length();
                //in case of deletions in UMI, we remove nucleotides from the last wildcard sequence
                //for the last UMI sequence, we take the whole sequence (in case of insertions)

                if( (lastWildcardEnd + currentWildcardLength) > oldWildcardMappedBarcode.length() || wildCardToFill == 1) 
                {
                    currentWildcardLength = oldWildcardMappedBarcode.length() - lastWildcardEnd;
                }

                std::string wildCardString = oldWildcardMappedBarcode.substr(lastWildcardEnd, currentWildcardLength);

                barcodeList.push_back(wildCardString);
                wildCardToFill -= 1;

                lastWildcardEnd += currentWildcardLength; // add the lengths of barcodes to the next offset position
            }
        }
        //add this match to the BarcodeMapping
        //barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        barcodeList.push_back(barcode);
    }
    //if the last barcode was a WIldcard that still has to be added:
    //BE CAREFUL: for now this means there can be ONLY ONE UMI at the end of a sequence
    if(wildCardToFill)
    {
        assert(wildCardToFill <= 1);
        wildCardToFill = 0; // unnecessary, still left to explicitely set to false
        std::string oldWildcardMappedBarcode = seq.first.substr(old_offset, seq.first.length() - old_offset);
        //barcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
        barcodeList.push_back(oldWildcardMappedBarcode);
    }

    barcodeMap.addVector(barcodeList);

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

bool MapEachBarcodeSequentiallyPolicyPairwise::map_forward(const std::string& seq, const input& input, 
                                                           BarcodePatternVectorPtr barcodePatterns,
                                                           fastqStats& stats,
                                                           std::vector<std::string>& barcodeList,
                                                           uint& barcodePosition,
                                                           int& score_sum)
{

    //iterate over BarcodeMappingVector
    int offset = 0;

    int old_offset = offset;
    unsigned int wildCardToFill = 0;
    int wildCardLength = 0, differenceInBarcodeLength = 0;
    for(BarcodePatternVector::iterator patternItr = barcodePatterns->begin(); 
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
            wildCardLength = (*patternItr)->get_patterns().at(0).length();
            offset += wildCardLength;
            wildCardToFill += 1;
            continue;
        }

        //for every barcodeMapping element find a match
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
        int start=0, end=0, score = 0;
        bool startCorrection = false;
        if(wildCardToFill){startCorrection = true;} // sart correction checks if we have to move our mapping window to the 5' direction
        
        // in case the sequence ends in the middle of the next barcode
        bool seqToShort = check_if_seq_too_short(offset, seq);
        if(seqToShort)
        {
            return true;
        }

        //if we did not match a pattern
        if(!(*patternItr)->match_pattern(seq, offset, start, end, score, barcode, differenceInBarcodeLength, startCorrection, false))
        {
            return false;
        }

        //set the length difference after barcode mapping
        //for the case of insertions inside the barcode sequence set the difference explicitely to zero 
        //(we only focus on deletions that we can not distinguish from substitutions)
        //differenceInBarcodeLength = barcode.length() - (end);
        if(differenceInBarcodeLength<0){differenceInBarcodeLength=0;}
        std::string sequence = seq.substr(offset + start, end-start);
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
            std::string oldWildcardMappedBarcode = seq.substr(old_offset, (offset+start-end) - old_offset);
            //barcodeMap.emplace_back(uniqueChars.  (oldWildcardMappedBarcode));
            unsigned int lastWildcardEnd = 0; //in case we have several wildcards and need to know old offset
            // to subset new string
            while(wildCardToFill != 0)
            {
                BarcodePatternVector::iterator wildCardIt = patternItr;
                int pos = -1 * wildCardToFill;
                std::advance(wildCardIt, pos);
                int currentWildcardLength = (*wildCardIt)->get_patterns().at(0).length();
                //in case of deletions in UMI, we remove nucleotides from the last wildcard sequence
                //for the last UMI sequence, we take the whole sequence (in case of insertions)
                if( (lastWildcardEnd + currentWildcardLength) > oldWildcardMappedBarcode.length() || wildCardToFill == 1) 
                {
                    currentWildcardLength = oldWildcardMappedBarcode.length() - lastWildcardEnd;
                }

                std::string wildCardString = oldWildcardMappedBarcode.substr(lastWildcardEnd, currentWildcardLength);
                barcodeList.push_back(wildCardString);
                wildCardToFill -= 1;

                lastWildcardEnd += currentWildcardLength; // add the lengths of barcodes to the next offset position
                ++barcodePosition; //increase the count of found positions
            }
        }

        //add this match to the BarcodeMapping
        //barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        barcodeList.push_back(barcode);
        ++barcodePosition; //increase the count of found positions
    }
    //if the last barcode was a WIldcard that still has to be added
    if(wildCardToFill)
    {
        wildCardToFill = 0; // unnecessary, still left to explicitely set to false
        std::string oldWildcardMappedBarcode = seq.substr(old_offset, seq.length() - old_offset);
        //barcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
        barcodeList.push_back(oldWildcardMappedBarcode);
        ++barcodePosition; //increase the count of found positions
    }

    return true;
}

bool MapEachBarcodeSequentiallyPolicyPairwise::map_reverse(const std::string& seq, const input& input, 
                                                           BarcodePatternVectorPtr barcodePatterns,
                                                           fastqStats& stats,
                                                           std::vector<std::string>& barcodeList,
                                                           uint& barcodePosition,
                                                           int& score_sum)
{
    //iterate over BarcodeMappingVector
    int offset = 0;

    int old_offset = offset;
    unsigned int wildCardToFill = 0;
    int wildCardLength = 0, differenceInBarcodeLength = 0;
    //iterate reverse through patterns
    for(BarcodePatternVector::reverse_iterator patternItr = barcodePatterns->rbegin(); 
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
            wildCardLength = (*patternItr)->get_patterns().at(0).length();
            offset += wildCardLength;
            wildCardToFill += 1;
            continue;
        }

        //for every barcodeMapping element find a match
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
        int start=0, end=0, score = 0;
        bool startCorrection = false;
        if(wildCardToFill){startCorrection = true;} // sart correction checks if we have to move our mapping window to the 5' direction
        
        // in case the sequence ends in the middle of the next barcode
        bool seqToShort = check_if_seq_too_short(offset, seq);
        if(seqToShort)
        {
            return true;
        }

        //map each pattern with reverse complement
        if(!(*patternItr)->match_pattern(seq, offset, start, end, score, barcode, differenceInBarcodeLength, startCorrection, true))
        {
            return false;
        }

        //set the length difference after barcode mapping
        //for the case of insertions inside the barcode sequence set the difference explicitely to zero 
        //(we only focus on deletions that we can not distinguish from substitutions)
        //we only substract the end, bcs we only consider barcode length differences at the end of the barcode
        //differenceInBarcodeLength = barcode.length() - (end);
        if(differenceInBarcodeLength<0){differenceInBarcodeLength=0;}
        std::string sequence = seq.substr(offset + start, end-start);
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
            std::string oldWildcardMappedBarcode = seq.substr(old_offset, (offset+start-end) - old_offset);

            //barcodeMap.emplace_back(uniqueChars.  (oldWildcardMappedBarcode));
            unsigned int lastWildcardEnd = 0; //in case we have several wildcards and need to know old offset
            // to subset new string
            while(wildCardToFill != 0)
            {
                BarcodePatternVector::reverse_iterator wildCardIt = patternItr;
                int pos = -1 * wildCardToFill;
                std::advance(wildCardIt, pos);
                int currentWildcardLength = (*wildCardIt)->get_patterns().at(0).length();
                //in case of deletions in UMI, we remove nucleotides from the last wildcard sequence
                //for the last UMI sequence, we take the whole sequence (in case of insertions)
                if( (lastWildcardEnd + currentWildcardLength) > oldWildcardMappedBarcode.length() || wildCardToFill == 1) 
                {
                    currentWildcardLength = oldWildcardMappedBarcode.length() - lastWildcardEnd;
                }

                std::string wildCardString = oldWildcardMappedBarcode.substr(lastWildcardEnd, currentWildcardLength);

                //in reverse mapping we have to make reverse complement of wildcard sequence
                std::string reverseComplimentbarcode = Barcode::generate_reverse_complement(wildCardString);

                barcodeList.push_back(reverseComplimentbarcode);
                wildCardToFill -= 1;

                lastWildcardEnd += currentWildcardLength; // add the lengths of barcodes to the next offset position
                ++barcodePosition; //increase the count of found positions
            }
        }



        //add this match to the BarcodeMapping
        //barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        barcodeList.push_back(barcode);
        ++barcodePosition; //increase the count of found positions
    }
    //if the last barcode was a WIldcard that still has to be added
    if(wildCardToFill)
    {
        wildCardToFill = 0; // unnecessary, still left to explicitely set to false
        std::string oldWildcardMappedBarcode = seq.substr(old_offset, seq.length() - old_offset);
        //barcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
        barcodeList.push_back(oldWildcardMappedBarcode);
        ++barcodePosition; //increase the count of found positions
    }

    return true;
}

bool MapEachBarcodeSequentiallyPolicyPairwise::combine_mapping(DemultiplexedReads& barcodeMap,
                                                               const BarcodePatternVectorPtr& barcodePatterns,
                                                               std::vector<std::string>& barcodeListFw,
                                                               const uint& barcodePositionFw,
                                                               const std::vector<std::string>& barcodeListRv,
                                                               const uint& barcodePositionRv,
                                                               fastqStats& stats,
                                                               int& score_sum,
                                                               std::pair<const std::string&, const std::string&> seq)
{
    //check that we span the whole sequence (except constant regions)
    int patternNum = barcodePatterns->size();

    //if positions are next to each other just return
    if(patternNum == (barcodePositionFw + barcodePositionRv))
    {
        for(std::vector<std::string>::const_reverse_iterator rvBarcodeIt = barcodeListRv.crbegin(); rvBarcodeIt != barcodeListRv.crend(); ++rvBarcodeIt)
        {
            barcodeListFw.push_back({*rvBarcodeIt});
        }
    }
    //if they r overlapping, check all overlapping ones are the same
    if(patternNum < (barcodePositionFw + barcodePositionRv))
    {
        //assert all overlapping barcodes r the same
        int start = (patternNum - barcodePositionRv); // this is the idx of the first shared barcode
        int end = barcodePositionFw - 1; //this is the idnex of the last shared barcode
        int j = barcodeListRv.size() - 1;
        for(int i = start; i <= end; ++i)
        {
            if(barcodeListFw.at(i) != barcodeListRv.at(j))
            {
                ++stats.noMatches;
                return false;
            }
            --j;
        }

        int overlap = (end - start) + 1;
        int missingPatternStartIdx = (barcodeListRv.size() - overlap) -1;
        //combine them to one vector
        for(int i = missingPatternStartIdx; i >= 0; --i)
        {
            barcodeListFw.push_back(barcodeListRv.at(i));
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
            if(!(barcodePatterns->at(i)->is_constant()))
            {
                ++stats.noMatches;
                return false;
            }
            barcodeListFw.push_back(barcodePatterns->at(i)->get_patterns().at(0));
        }

        //add reverse patterns
        for(std::vector<std::string>::const_reverse_iterator rvBarcodeIt = barcodeListRv.crbegin(); rvBarcodeIt != barcodeListRv.crend(); ++rvBarcodeIt)
        {
            barcodeListFw.push_back({*rvBarcodeIt});
        }
    }

    //add the final barcodeList
    barcodeMap.addVector(barcodeListFw);

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

bool MapEachBarcodeSequentiallyPolicyPairwise::split_line_into_barcode_patterns(std::pair<const std::string&, const std::string&> seq, const input& input, 
                                        DemultiplexedReads& barcodeMap, 
                                        BarcodePatternVectorPtr barcodePatterns,
                                        fastqStats& stats)
{

    int score_sum = 0;
    //map forward reads barcodeList contains the stored barcodes, 
    //barcodePosition is the psoiton of the last mapped barcode
    std::vector<std::string> barcodeListFw;
    uint barcodePositionFw = 0;
    bool fwBool = map_forward(seq.first, input, barcodePatterns, stats, barcodeListFw, barcodePositionFw, score_sum);

    std::vector<std::string> barcodeListRv;
    uint barcodePositionRv = 0;
    bool rvBool = map_reverse(seq.second, input, barcodePatterns, stats, barcodeListRv, barcodePositionRv, score_sum);

    combine_mapping(barcodeMap, barcodePatterns, barcodeListFw, barcodePositionFw, barcodeListRv, barcodePositionRv, stats, score_sum, seq);

    return true;
}

template <typename MappingPolicy, typename FilePolicy>
bool Mapping<MappingPolicy, FilePolicy>::demultiplex_read(std::pair<const std::string&, const std::string&> seq, const input& input, 
                                                          std::atomic<unsigned long long>& count, const unsigned long long& totalReadCount,
                                                          bool guideMapping)
{
    //split line into patterns (barcodeMap, barcodePatters, stats are passed as reference or ptr)
    //and can be read by each thread, "addValue" method for barcodeMap is thread safe also for concurrent writing
    bool result;
    if(!guideMapping)
    {
        result = this->split_line_into_barcode_patterns(seq, input, barcodeMap, barcodePatterns, stats);
    }
    else
    {
        result = this->split_line_into_barcode_patterns(seq, input, guideBarcodeMap, guideBarcodePatterns, stats);
        --stats.noMatches;
    }

    //update status bar
    ++count;
    if(count%1000==0 && totalReadCount!=ULLONG_MAX) //update at every 1,000th entry
    {
        double perc = count/(double)totalReadCount;
        std::lock_guard<std::mutex> guard(*printProgressLock);
        printProgress(perc);
    }

    return(result);
}

void MapAroundConstantBarcodesAsAnchorPolicy::map_pattern_between_linker(const std::string& seq, const int& oldEnd, 
                                                                         const int& start,
                                                                         BarcodePatternVectorPtr barcodePatterns,
                                                                         std::vector<std::string>& barcodeList,
                                                                         int& barcodePosition, int& skippedBarcodes)
{
    //we have to map all barcodes that we missed
    int skipend=oldEnd;
    std::string skippedBarcodeString = seq.substr(skipend, start-skipend);
    int skipStringNewOffset = 0;
    for(int i = barcodePosition-skippedBarcodes; i < barcodePosition ; ++i)
    {
        BarcodePatternVector::iterator skippedPatternItr = barcodePatterns->begin() + i;
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

bool MapAroundConstantBarcodesAsAnchorPolicy::split_line_into_barcode_patterns(std::pair<const std::string&, const std::string&> seq, const input& input, 
                                        DemultiplexedReads& barcodeMap, 
                                        BarcodePatternVectorPtr barcodePatterns,
                                        fastqStats& stats)
{
    int offset = 0;
    int oldEnd = 0;

    std::vector<std::string> barcodeList;

    //firstly map each constant barcode
    int barcodePosition = 0;
    int skippedBarcodes = 0;
    for(BarcodePatternVector::iterator patternItr = barcodePatterns->begin(); 
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
            
        std::string subStringToSearchBarcodes = seq.first.substr(oldEnd, seq.first.length() - oldEnd);

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
            map_pattern_between_linker(seq.first,oldEnd, start, barcodePatterns, barcodeList, barcodePosition, skippedBarcodes);
            //std::string skippedBarcodeString = seq.first.substr(oldEnd, start-oldEnd);
            //barcodeList.push_back(skippedBarcodeString);
        }

        //2.) write the newly mapped constant barcode
        barcodeList.push_back(barcode);

        //set paramters after constant mapping
        oldEnd += end;
        ++barcodePosition;
        skippedBarcodes = 0;
    }

    if(skippedBarcodes > 0)
    {
        int start = seq.first.length();
        map_pattern_between_linker(seq.first,oldEnd, start, barcodePatterns, barcodeList, barcodePosition, skippedBarcodes);
        //std::string skippedBarcodeString = seq.first.substr(oldEnd, start-oldEnd);
        //barcodeList.push_back(skippedBarcodeString);
    }

    barcodeMap.addVector(barcodeList);
    ++stats.perfectMatches;

    return true;
}

template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::run_mapping(const input& input)
{
    std::cout << "START DEMULTIPLEXING\n";

    //generate a pool of threads
    boost::asio::thread_pool pool(input.threads); //create thread pool

    //read line by line and add to thread pool
    FilePolicy::init_file(input.inFile, input.reverseFile);
    std::pair<std::string, std::string> line;
    std::atomic<unsigned long long> lineCount = 0; //using atomic<int> as thread safe read count
    unsigned long long totalReadCount = FilePolicy::get_read_number();

    while(FilePolicy::get_next_line(line))
    {
        //handing over only lineCount as reference, everything else will be copied (Mapping object as handed overr as this-pointer)
        //be aware: in default function do not handle guide reads, this is part of the overwritten function in Demultiplexing tool
        boost::asio::post(pool, std::bind(&Mapping::demultiplex_read, this, line, input, std::ref(lineCount), totalReadCount, false));
    }
    pool.join();
    printProgress(1); std::cout << "\n"; // end the progress bar
    if(totalReadCount != ULLONG_MAX)
    {
        std::cout << "=>\t READS WITH A MATCHED BARCODE: " << std::to_string((unsigned long long)(100*(stats.perfectMatches)/(double)totalReadCount)) 
                << "% | MODERATE MATCHES: " << std::to_string((unsigned long long)(100*(stats.moderateMatches)/(double)totalReadCount))
                << "% | Linker sequences mapped non sequentially (e.g. same linker sequences): " << std::to_string((unsigned long long)(100*(stats.noMatches)/(double)totalReadCount)) << "%\n";
    }
    FilePolicy::close_file();
}


template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::run(const input& input)
{
    //generate barcode objects that store all the information(e.g. allowed barcodes, mismatches)
    generate_barcode_patterns(input);

    //
    if(input.writeStats)
    {
        initializeStats();
    }

    //run over all reads and sequentially map every barcode
    run_mapping(input);
}

template class Mapping<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy>;
template class Mapping<MapEachBarcodeSequentiallyPolicyPairwise, ExtractLinesFromFastqFilePolicyPairedEnd>;
template class Mapping<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy>;
template class Mapping<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromTxtFilesPolicy>;
template class Mapping<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromFastqFilePolicy>;
