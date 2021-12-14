#include <iostream>

#include "mapping.hpp"
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

namespace
{

void write_failed_lines(const input& input, const std::vector<std::string>& failedLines)
{
    std::string output = input.outFile;
    std::ofstream outputFile;

    //write real sequences that map to barcodes
    std::size_t found = output.find_last_of("/");
    std::string outputFail;
    if(found == std::string::npos)
    {
        outputFail = "FailedLines_" + output;
    }
    else
    {
        outputFail = output.substr(0,found) + "/" + "FailedLines_" + output.substr(found+1);
    }
    outputFile.open (outputFail, std::ofstream::app);
    for(int i = 0; i < failedLines.size(); ++i)
    { 
        outputFile << failedLines.at(i) << "\n";
    }
    outputFile.close();
}

void writeStats(std::string output, const fastqStats& stats)
{
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");
    if(found == std::string::npos)
    {
        output = "StatsBarcodeMappingErrors_" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "StatsBarcodeMappingErrors_" + output.substr(found+1);
    }
    std::remove(output.c_str());
    outputFile.open (output, std::ofstream::app);

    for(std::pair<std::string, std::vector<int> > mismatchDictEntry : stats.mapping_dict)
    {
        outputFile << mismatchDictEntry.first << "\t";
        for(int mismatchCountIdx = 0; mismatchCountIdx < mismatchDictEntry.second.size(); ++mismatchCountIdx)
        {
            outputFile << mismatchDictEntry.second.at(mismatchCountIdx);
            if(mismatchCountIdx!=mismatchDictEntry.second.size()-1){outputFile << "\t";}
        }
        outputFile << "\n";
    }

    outputFile.close();
}

/*void writeCurrentResults(BarcodeMappingVector& barcodes, BarcodeMappingVector& realBarcodes, const input& input, 
                         std::vector<std::string>& failedLines)
{
    write_file(input, barcodes, realBarcodes);
    write_failed_lines(input, failedLines);

    barcodes.clear();
    realBarcodes.clear();
    failedLines.clear();
}*/

}

bool MapEachBarcodeSequentiallyPolicy::check_if_seq_too_short(const int& offset, const std::string& seq)
{
    if(offset >= seq.length())
    {
        return true;
    }

    return false;
}

template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::parse_barcode_data(const input& input, std::vector<std::pair<std::string, char> >& patterns, std::vector<int>& mismatches, std::vector<std::vector<std::string> >& varyingBarcodes)
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
                            std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
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
            exit(1);
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
    int variableBarcodeIdx = 0; // index for variable barcode sequences is shorter than the vector of patterns in total
    for(int i=0; i < patterns.size(); ++i)
    {
        if(patterns.at(i).second=='v')
        {
            VariableBarcode barcode(varyingBarcodes.at(variableBarcodeIdx), mismatches.at(i));
            std::shared_ptr<VariableBarcode> barcodePtr(std::make_shared<VariableBarcode>(barcode));
            barcodeVector.emplace_back(barcodePtr);
            ++variableBarcodeIdx;
        }
        else if(patterns.at(i).second=='c')
        {
            ConstantBarcode barcode(patterns.at(i).first, mismatches.at(i));
            std::shared_ptr<ConstantBarcode> barcodePtr(std::make_shared<ConstantBarcode>(barcode));
            barcodeVector.emplace_back(barcodePtr);
        }
        else if(patterns.at(i).second=='w')
        {
            WildcardBarcode barcode(patterns.at(i).first, mismatches.at(i));
            std::shared_ptr<WildcardBarcode> barcodePtr(std::make_shared<WildcardBarcode>(barcode));
            barcodeVector.emplace_back(barcodePtr);
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
    return patterns;
}

//END INITIALIZATION HELPER FUNCTIONS















//apply all BarcodePatterns to a single line to generate a vector strings (the Barcodemapping)
//realBarcodeMap contains the actual string in the seauence that gets a barcode assigned
bool MapEachBarcodeSequentiallyPolicy::split_line_into_barcode_patterns(const std::string& seq, const input& input, 
                                        DemultiplexedReads& barcodeMap, 
                                        BarcodePatternVectorPtr barcodePatterns)
{
    MappingResult result;
    std::vector<std::string> barcodeList;

    //iterate over BarcodeMappingVector
    int offset = 0;
    int score_sum = 0;

    int old_offset = offset;
    bool wildCardToFill = false;
    int wildCardLength = 0, differenceInBarcodeLength = 0;
    for(BarcodePatternVector::iterator patternItr = barcodePatterns->begin(); 
        patternItr < barcodePatterns->end(); 
        ++patternItr)
    {
        //if we have a wildcard skip this matching, we match again the next sequence
        if((*patternItr)->is_wildcard())
        {
            old_offset = offset;
            wildCardLength = (*patternItr)->get_patterns().at(0).length();
            offset += wildCardLength;
            wildCardToFill = true;
            continue;
        }
        //for every barcodeMapping element find a match
        std::string barcode = ""; //the actual real barcode that we find (mismatch corrected)
        int start=0, end=0, score = 0;
        bool startCorrection = false;
        if(wildCardToFill){startCorrection = true;} // sart correction checks if we have to move our mapping window to the 5' direction
        // could happen in the case of deletions in the UMI sequence...
        bool seqToShort = check_if_seq_too_short(offset, seq);
        if(seqToShort && input.analyseUnmappedPatterns)
        {
            break;
        }
        else if(seqToShort)
        {
            return false;
        }
        if(!(*patternItr)->match_pattern(seq, offset, start, end, score, barcode, differenceInBarcodeLength, startCorrection))
        {
            //++stats.noMatches;
            if(!input.analyseUnmappedPatterns)
            {
                return false;
            }
            else
            {
                barcode = "";
                start = 0;
                end = 0;
            }
        }
        //set the length difference after barcode mapping
        //for the case of insertions inside the barcode sequence set the difference explicitely to zero 
        //(we only focus on deletions that we can not distinguish from substitutions)
        differenceInBarcodeLength = barcode.length() - (end-start);
        if(differenceInBarcodeLength<0){differenceInBarcodeLength=0;}
        std::string sequence = seq.substr(offset + start, end-start);
        offset += end;
        score_sum += score;

        if(!input.analyseUnmappedPatterns){assert(barcode != "");}
        //add barcode data to statistics dictionary
        int dictvectorIndex = ( (score <= (*patternItr)->mismatches) ? (score) : ( ((*patternItr)->mismatches) + 1) );
        
        //if(!input.analyseUnmappedPatterns)
        //{
         //   ++stats.mapping_dict[realBarcode].at(dictvectorIndex);
        //}
        
        //squeeze in the last wildcard match if there was one 
        //(we needed both matches of neighboring barcodes to define wildcard boundaries)
        if(wildCardToFill == true)
        {
            wildCardToFill = false;
            startCorrection = false;
            std::string oldWildcardMappedBarcode = seq.substr(old_offset, (offset+start-end) - old_offset);
            //barcodeMap.emplace_back(uniqueChars.  (oldWildcardMappedBarcode));
            barcodeList.push_back(oldWildcardMappedBarcode);
            result.sequences.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
        }
        //add this match to the BarcodeMapping
        //barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        barcodeList.push_back(barcode);
        result.sequences.emplace_back(std::make_shared<std::string>(sequence));
    }
    //if the last barcode was a WIldcard that still has to be added
    if(wildCardToFill == true)
    {
        wildCardToFill = false; // unnecessary, still left to explicitely set to false
        std::string oldWildcardMappedBarcode = seq.substr(old_offset, seq.length() - old_offset);
        //barcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
        barcodeList.push_back(oldWildcardMappedBarcode);
        result.sequences.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
    }

    barcodeMap.addVector(barcodeList);
    result.score = score_sum;

    return true;
}

template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::demultiplex_read(const std::string& seq, const input& input, 
                                                          std::atomic<int>& count, const int& totalReadCount)
{

    // split line into patterns
    this->split_line_into_barcode_patterns(seq, input, 
                                     barcodeMap, 
                                      barcodePatterns);

    //update status bar
    ++count;
    if(count%1000 == 0) //update at every 1,000th entry
    {
        double perc = count/(double)totalReadCount;
        std::lock_guard<std::mutex> guard(*mappingLock);
        printProgress(perc);
    }
}

template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::run_mapping(const input& input)
{
    std::cout << "START DEMULTIPLEXING\n";

    //generate a pool of threads
    boost::asio::thread_pool pool(input.threads); //create thread pool

    //read line by line and add to thread pool
    FilePolicy::init_file(input.inFile);
    std::string line;
    std::atomic<int> lineCount = 0; //using atomic<int> as thread safe read count
    int totalReadCount = numberOfReads(input.inFile);

    while(FilePolicy::get_next_line(line))
    {
        //if we processed a number of lines wait for threads to finish
        boost::asio::post(pool, std::bind(&Mapping::demultiplex_read, this, line, input, std::ref(lineCount), totalReadCount));
    }
    pool.join();
    printProgress(1); std::cout << "\n"; // end the progress bar

    std::cout << barcodeMap.size() << "\n";
}


template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::run(const input& input)
{
    generate_barcode_patterns(input);

    run_mapping(input);
}



















/*
bool MapAroundConstantBarcodesAsAnchorPolicy::split_line_into_barcode_patterns(const std::string& seq, const input& input, BarcodeMapping& barcodeMap,
                                      BarcodePatternVectorPtr barcodePatterns)
{
    int offset = 0;

    int oldEnd = 0;

    //map each constant barcode
    int barcodePosition = 0;
    for(BarcodePatternVector::iterator patternItr = barcodePatterns->begin(); 
        patternItr < barcodePatterns->end(); 
        ++patternItr)
    {
        //exclude non constant
        if(!(*patternItr)->is_constant()){++barcodePosition;continue;}
        
        //map next constant region
        int start=0, end=0, score = 0, differenceInBarcodeLength = 0;
        std::string realBarcode = ""; //the actual real barcode that we find (mismatch corrected)
        if(!(*patternItr)->match_pattern(seq, offset, start, end, score, realBarcode, differenceInBarcodeLength, false, true))
        {
            ++barcodePosition;
            continue;
        }

        std::string mappedBarcode = seq.substr(offset + start, end-start);

        if(start < oldEnd){return false;} //we have in this case barcodes that are non sequential

        int currentPosition = barcodeMap.size();
        for(int i = currentPosition; i <= barcodePosition; ++i)
        {
            std::string skippedBarcodes = seq.substr(oldEnd, start-oldEnd);
            //the first barcode which is not the constant barcode
            if(i ==  currentPosition)
            {
                barcodeMap.emplace_back(std::make_shared<std::string>(skippedBarcodes));
                //realBarcodeMap.emplace_back(std::make_shared<std::string>(skippedBarcodes));
            }
            //if it is the constant barcode
            else if(i == barcodePosition)
            {
                barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
                //realBarcodeMap.emplace_back(std::make_shared<std::string>(realBarcode));
            }
            //every other barcode
            else
            {
                barcodeMap.emplace_back(std::make_shared<std::string>(""));
                //realBarcodeMap.emplace_back(std::make_shared<std::string>(""));
            }

        }
        oldEnd = end;
        ++barcodePosition;
    }
    return true;
}

template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::map_pattern_to_fastq_lines(std::vector<std::string>& fastqLines, const input& input, BarcodeMappingVector& barcodes, 
                                BarcodeMappingVector& realBarcodes, fastqStats& stats, BarcodePatternVectorPtr barcodePatterns,
                                std::vector<std::string>& failedLines)
{
    for(const std::string& line : fastqLines)
    {
        BarcodeMapping barcode;
        if(MappingPolicy::split_line_into_barcode_patterns(line, input, barcode, barcodePatterns))
        {
            barcodes.push_back(barcode);
            //realBarcodes.push_back(realBarcode);
            barcode.clear();
        }
        else
        {
            failedLines.push_back(line);
        }
    }
}

template <typename MappingPolicy, typename FilePolicy>
void Mapping<MappingPolicy, FilePolicy>::ditribute_jobs_to_threads(const input& input, std::vector<std::string>& fastqLines)
{
    //split input lines into thread buckets
    std::vector<std::vector<std::string> > fastqLinesVector; //vector holding all the buckets of fastq lines to be analyzed by each thread
                                                            // each vector element is one bucket for one thread
    int element_number = fastqLines.size() / input.threads;
    std::vector<std::string>::iterator begin = fastqLines.begin();
    for(int i = 0; i < input.threads; ++i)
    {
        std::vector<std::string>::iterator end = ( i==input.threads-1 ? fastqLines.end() : begin + element_number);

        std::vector<std::string> tmpFastqLines(begin, end);
        fastqLinesVector.push_back(tmpFastqLines);
        begin = end;
    }
    //for every bucket call a thread
    //tmp variables to store thread results
    std::vector<std::string> failedLines;
    std::vector<std::vector<std::string>> failedLinesThreadlist(input.threads);
    std::vector<BarcodeMappingVector> barcodesThreadList(input.threads); // vector of mappedSequences to be filled
    std::vector<BarcodeMappingVector> realBarcodesThreadList(input.threads); // vector of mappedSequences to be filled with corrected matches
    std::vector<fastqStats> statsThreadList(input.threads, *fastqStatsPtr); // vector of stats to be filled
    std::vector<std::thread> workers;
    for (int i = 0; i < input.threads; ++i) {
        workers.push_back(std::thread(&Mapping<MappingPolicy, FilePolicy>::map_pattern_to_fastq_lines, this, std::ref(fastqLinesVector.at(i)), 
                        std::ref(input), std::ref(barcodesThreadList.at(i)), std::ref(realBarcodesThreadList.at(i)), 
                        std::ref(statsThreadList.at(i)), barcodePatterns, std::ref(failedLinesThreadlist.at(i))));

           // workers.push_back(std::thread(xxx, std::ref(fastqLinesVector.at(i)), std::ref(input)));
    }
    for (std::thread &t: workers) 
    {
        if (t.joinable()) {
            t.join();
        }
    }
    //combine thread data
    for (int i = 0; i < input.threads; ++i) 
    {
        //sequenceBarcodes.insert(sequenceBarcodes.end(), barcodesThreadList.at(i).begin(), barcodesThreadList.at(i).end());
        realBarcodes.insert(realBarcodes.end(), realBarcodesThreadList.at(i).begin(), realBarcodesThreadList.at(i).end());

        if(!input.analyseUnmappedPatterns)
        {
            fastqStatsPtr->perfectMatches += statsThreadList.at(i).perfectMatches;
            fastqStatsPtr->noMatches += statsThreadList.at(i).noMatches;
            fastqStatsPtr->moderateMatches += statsThreadList.at(i).moderateMatches;
            fastqStatsPtr->multiBarcodeMatch += statsThreadList.at(i).multiBarcodeMatch;
            failedLines.insert(failedLines.end(), failedLinesThreadlist.at(i).begin(), failedLinesThreadlist.at(i).end());
            //iterate for every thread over the dictionary of mismatches per barcode and fill the final stats data with it
            for(std::pair<std::string, std::vector<int> > mismatchDictEntry : statsThreadList.at(i).mapping_dict)
            {
                for(int mismatchCountIdx = 0; mismatchCountIdx < mismatchDictEntry.second.size(); ++mismatchCountIdx)
                {
                    fastqStatsPtr->mapping_dict[mismatchDictEntry.first].at(mismatchCountIdx) += statsThreadList.at(i).mapping_dict[mismatchDictEntry.first].at(mismatchCountIdx);
                }
            }
        }
    }
    //writeCurrentResults(sequenceBarcodes, realBarcodes, input, failedLines);
}*/

template class Mapping<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy>;
template class Mapping<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy>;
//template class Mapping<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromTxtFilesPolicy>;
