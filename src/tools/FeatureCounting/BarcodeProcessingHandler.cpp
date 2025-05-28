#include "BarcodeProcessingHandler.hpp"

#include <unordered_set>
#include <unordered_map>
#include <set>
#include <cstdlib>

double calcualtePercentages(std::vector<unsigned long long> groups, int num, double perc)
{
    double readCount = perc * num;
    int sumOfReads = 0;
    int sumBCs = 0;
    std::sort(groups.begin(), groups.end(), std::greater<unsigned long long>());
    for(auto el : groups)
    {
        sumOfReads += el;
        ++sumBCs;
        if(double(sumOfReads) >= readCount)
        {
            return((double)sumBCs/groups.size());
        }
    }
    assert(groups.size() != 0);
    return(-1);
}

bool endsWithTxt(const std::string& filename) 
{
    return filename.size() >= 4 && filename.substr(filename.size() - 4) == ".txt";
}

int isUMICol(const std::string& str) 
{
    std::regex pattern(R"(^(\d+)X$)");
    std::smatch match;
    
    if (std::regex_match(str, match, pattern)) 
    {
        return std::stoi(match[1]); 
    }
    return 0;
}

void parseVariableBarcodeFile(const std::string& file, std::unordered_map<int, std::vector<std::string>>& barcodeList, int colIdx)
{
    std::ifstream barcodeFileStream(file);

    std::string line;
    while (std::getline(barcodeFileStream, line)) 
    {
        std::string delimiter = ",";
        std::string seq;
        size_t pos = 0;
        std::vector<std::string> seqVector;
        while ((pos = line.find(delimiter)) != std::string::npos) 
        {
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
        barcodeList.insert(std::make_pair(colIdx, seqVector));
        seqVector.clear();
    }
}

std::vector<std::string> parseFeatureNames(const std::string& featureFile) {
    std::vector<std::string> result;
    std::stringstream ss(featureFile);
    std::string item;
    while (std::getline(ss, item, ',')) 
    {
        result.push_back(item);
    }
    return result;
}

void generateBarcodeDicts(const std::string& headerLine, const std::string& barcodeDir, std::string barcodeIndices, BarcodeInformation& barcodeIdData, 
                          std::vector<std::string>& proteinNamelist, bool parseAbBarcodes, const int& featureIdx, 
                          std::vector<std::string>* treatmentDict, const int& treatmentIdx, std::string umiIdx,
                          int umiMismatches)
{

    //parse UMI columns into temporary vector, when parsing header below check if the 
    //index that we parsed here is indeed a umI col (e.g., 10X)
    std::vector<int> tmpUmiIdices;
    if(umiIdx != "")
    {
        std::stringstream ssUmi;
        ssUmi.str(umiIdx);
        // parse all the indices that r used for the UMI (can be several)
        //when parsing the header make sure the umiIdices are indeed a random sequence
        while(ssUmi.good())
        {
            std::string umiSubstr;
            getline(ssUmi, umiSubstr, ',' );
            tmpUmiIdices.push_back(stoi(umiSubstr));
        }
    }

    //assign number of UMI mismatches
    barcodeIdData.umiMismatches = umiMismatches;

    //parse barcode file into a vector of a vector of all sequences
    std::vector<std::string> barcodeHeader;

    //iterate through the header of the input file
    std::stringstream ss(headerLine);
    std::string colName;
    //iterate over columns
    int colIdx = 0;

    //PROCESS HEADERS
    std::unordered_map<int, std::vector<std::string>> barcodeList; //maps column number -> all barcodes: e.g. 2 -> variables barcode sin BC1.txt, ...
    while (std::getline(ss, colName, '\t')) 
    {
        //STORE ALL HEADERS
        barcodeHeader.push_back(colName);

        //MAPPING OF COLIDX -> LIST OF ALL BARCODES (STRINGS)
        if(endsWithTxt(colName))
        {
            parseVariableBarcodeFile(barcodeDir + "/" + colName, barcodeList, colIdx);
        }
        else
        {
            //MAP OF COLIDX to UMI, UMILENGTH
            int umiLen = isUMICol(colName);
            bool isAssignedUmi;
            if(umiIdx == ""){isAssignedUmi = true;}
            else
            {
                isAssignedUmi = std::find(tmpUmiIdices.begin(), tmpUmiIdices.end(), colIdx) != tmpUmiIdices.end();
            }
            //the returned UMILEN is the length of the umi sequence, in case it is a UMI
            if(umiLen && isAssignedUmi)
            {
                barcodeIdData.umiIdx.push_back(colIdx);
                barcodeIdData.umiLength += umiLen;
            }
        }
        ++colIdx;
    }

    //if we have CombinatorialIndexing lines, parse the single cell IDs
    if(barcodeIndices != "")
    {
        //parse the indices of CI-barcodes (the index of the lines that store CI-barcodes)
        std::stringstream ss;
        ss.str(barcodeIndices);

        // barcodeIdData.scBarcodeIndices contains a list of all COLIDX for sc assignment
        while(ss.good())
        {
            std::string substr;
            getline(ss, substr, ',' );
            barcodeIdData.scBarcodeIndices.push_back(stoi(substr));
        }

        std::cout << "Assinging single-cells according to columns: ";
        for(int i = 0; i < barcodeIdData.scBarcodeIndices.size()-1; ++i)
        {
            std::cout << barcodeHeader.at(barcodeIdData.scBarcodeIndices.at(i)) << ",";
        }
        std::cout << barcodeHeader.at(barcodeIdData.scBarcodeIndices.back());
        std::cout << "\n";
        
        //make map colIdx -> (map: barcode string -> number)
        //iterate over colIdx for scID
        for(int scIdx : barcodeIdData.scBarcodeIndices)
        {
            int barcodeCount = 0;
            std::unordered_map<std::string, int> barcodeMap;
            //map all barcode strings to numbers
            for(const std::string& barcodeEntry : barcodeList.at(scIdx))
            {
                barcodeMap.insert(std::pair<std::string, int>(barcodeEntry, barcodeCount));
                ++barcodeCount;
            }
            barcodeIdData.barcodeIdMaps.push_back(barcodeMap);
        }
    }

    //print assigned grouping index
    if(treatmentIdx!=-1)
    {
        std::cout << "Assigning treatment to column: " << barcodeHeader.at(treatmentIdx) << "\n";
        barcodeIdData.treatmentIdx = treatmentIdx;    
    }

    //print the used feature index
    std::cout << "Assigning feature to column: " << barcodeHeader.at(featureIdx) << "\n";
    barcodeIdData.featureIdx = featureIdx;

    //print the used UMI indices
    std::cout << "Using following columns as UMIs: ";
    for(int colIdxForUMI : barcodeIdData.umiIdx)
    {
        std::cout << barcodeHeader.at(colIdxForUMI) << ",";
    }
    std::cout << "\n";
    std::cout << "Aligning all UMI-sequences with " << std::to_string(barcodeIdData.umiMismatches) << " mismatches.\n";

    //assign the final dictionaries of variable barcodes
    if(parseAbBarcodes)
    {
        proteinNamelist = barcodeList.at(featureIdx);
    }
    if(treatmentIdx!=-1)
    {
        *treatmentDict = barcodeList.at(treatmentIdx);
    }

}

void BarcodeProcessingHandler::parse_barcode_file(const std::string fileName, const int& thread)
{
    unsigned long long totalReads = totalNumberOfLines(fileName);
    unsigned long long currentReads = 0;
    //open gz file
    if(!endWith(fileName,".gz"))
    {
        std::cerr << "Input file must be gzip compressed\n";
        exit(EXIT_FAILURE);
    }
    std::ifstream file(fileName, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(file);
    std::istream instream(&inbuf);
    
    parseBarcodeLines(&instream, totalReads, currentReads);

    file.close();
}

void BarcodeProcessingHandler::parseBarcodeLines(std::istream* instream, const unsigned long long& totalReads, unsigned long long& currentReads)
{
    std::string line;
    std::cout << "STEP[1/3]\t(READING ALL LINES INTO MEMORY)\n";
    int elements = 0; //check that each row has the correct number of barcodes
    unsigned long long readCount = 0;
    while(std::getline(*instream, line))
    {
        //Skip the header line
        if(currentReads==0)
        {
            std::stringstream ss(line);
            std::string item;
            while (std::getline(ss, item, '\t')) 
            {
                elements++;
            }
            ++currentReads; 
            continue;
        }

        add_line_to_temporary_data(line, elements, readCount);   

        double perc = currentReads/ (double)totalReads;
        ++currentReads;
        printProgress(perc);        
    }

    result.set_total_reads(currentReads-1); //minus header line
    result.set_total_ab_reads(readCount);

    printProgress(1);
    std::cout << "\n";
}

void BarcodeProcessingHandler::add_line_to_temporary_data(const std::string& line, const int& elements, unsigned long long& readCount)
{
    //split the line into barcodes
    std::vector<std::string> result;
    std::stringstream ss;
    ss.str(line);
    std::string substr;

    while(getline( ss, substr, '\t' ))
    {
        if(substr != ""){result.push_back( substr );}
    }
    if(result.size() != elements)
    {
        std::cout << "WARNING in barcode file, following row has not the correct number of sequences: " << line << "\n";
        return;
    }

    //hand over the UMI string, ab string, singleCellstring (concatenation of CIbarcodes)
    std::vector<std::string> ciBarcodes;
    //CiBarcodes are added in order of the scBarcodeIndices, this order must always be respected!!!
    for(int i : barcodeInformation.scBarcodeIndices)
    {
        ciBarcodes.push_back(result.at(i));
    }
    std::string singleCellIdx = generateSingleCellIndexFromBarcodes(ciBarcodes);

    std::string featureName = "";
    featureName = rawData.getFeatureName(result.at(barcodeInformation.featureIdx));
    
    std::string treatment = "";
    if(barcodeInformation.treatmentIdx != -1)
    {
        treatment = rawData.getTreatmentName(result.at(barcodeInformation.treatmentIdx));
    }

    ++readCount;
    const char* umiSeq;
    //if there is a UMI and also we should filter reads by the fact that a UMI should belong only to one SC-AB
    //the also create a UMI-SCAB Dict for filtering
    //(this is only useful if we expected the data to be extremely noisy or so shallow that there no
    //UMI-clashes: e.g. for debugging of CI experiments with many barcode recombinations to reduce erroneous reads)
    if(!barcodeInformation.umiIdx.empty() && umiRemoval)
    {
        std::string umiSeqString;
        for(int idx : barcodeInformation.umiIdx)
        {
            std::string tmpUmi = result.at(idx).c_str();
            umiSeqString = umiSeqString + tmpUmi;
        }
        umiSeq = umiSeqString.c_str();

        rawData.add_to_umiDict(umiSeq, featureName, singleCellIdx, treatment);
    }
    //otherwise add reads directly to dict of ScAb to reads
    else
    {
        rawData.add_to_scAbDict("", featureName, singleCellIdx, treatment);
    }
}

std::string BarcodeProcessingHandler::generateSingleCellIndexFromBarcodes(const std::vector<std::string>& ciBarcodes)
{
    std::string scIdx;
    if(scIdString)
    {
        for(int i=0; i < ciBarcodes.size(); ++i)
        {
            scIdx += ciBarcodes.at(i);
        }
        return scIdx;
    }
    for(int i = 0; i < barcodeInformation.scBarcodeIndices.size(); ++i)
    {
        //first is the column idx of the barcode, second is the actual barcode that we want to map to a number
        int tmpIdx = (barcodeInformation.barcodeIdMaps.at(i)).at(ciBarcodes.at(i));
        scIdx += std::to_string(tmpIdx);
        if(i < ciBarcodes.size() - 1)
        {
            scIdx += ".";
        }
    }

    return scIdx;
}

//all reads of the same UMI are combined -> written to a dict for an AB of a unique cell (this read is stored only once, there can already be seen
//as a UMI collapsing step)
void BarcodeProcessingHandler::markReadsWithNoUniqueUmi(const std::vector<umiDataLinePtr>& uniqueUmis,
                                                        std::atomic<unsigned long long>& count,
                                                        const unsigned long long& totalCount)
{
    //count how often we see which AB-SC combinations for this certain UMI
    std::unordered_map<std::string, unsigned long long> umiCountMap;
    std::unordered_map<std::string, umiDataLinePtr> umiReadMap; //mapping the unique ID to the first occuring actual read

    unsigned long long totalReadCount = uniqueUmis.size();
    for(unsigned long long i = 0; i < uniqueUmis.size(); ++i)
    {
        std::string uniqueID = std::string(uniqueUmis.at(i)->scID) + std::string(uniqueUmis.at(i)->abName);
        std::unordered_map<std::string, unsigned long long>::iterator umiCountMapIt = umiCountMap.find(uniqueID);
        if( umiCountMapIt != umiCountMap.end())
        {
            ++(umiCountMapIt->second);
        }
        else
        {
            umiCountMap.insert(std::make_pair(uniqueID, 1));
            umiReadMap.insert(std::make_pair(uniqueID, uniqueUmis.at(i)));
        }
    }
    //check there is a single cell + AB combination representing more than 90% of the UMI reads
    std::vector<std::string> realSingleCellIDVec;
    bool realSingleCellExists = false;
    for(auto singleCellCountPair : umiCountMap)
    {
        double singleCellPerc = (double)singleCellCountPair.second/totalReadCount;
        if( singleCellPerc >= umiFilterThreshold) //default = 0.9
        {
            realSingleCellIDVec.push_back(singleCellCountPair.first);
            realSingleCellExists = true;
            if(umiFilterThreshold>0.5){break;} //in this case we have only one entry in the vector
        }
    }
    //Collapse UMIs, remove false reads, map a single cell class name to single cells
    unsigned long long readsWithNoClass = 0;
    unsigned long long readsToKeep = 0;
    if(realSingleCellExists)
    {
        //delete only the <=10% 'false' reads
        for(auto it = umiReadMap.begin(); it != umiReadMap.end(); it++)
        {
            std::string uniqueID = it->first;
            if(std::find(realSingleCellIDVec.begin(), realSingleCellIDVec.end(), uniqueID) != realSingleCellIDVec.end())//add to new DIct
            {
                    const char* className = nullptr;
                    if(rawData.check_class())
                    {
                        className = rawData.get_sc_class_name(it->second->scID);
                        if(className == nullptr && scMustHaveClass)
                        {
                            //single cell has no class name
                            readsWithNoClass += umiCountMap.at(uniqueID);
                            if(umiFilterThreshold>0.5){break;}
                            else{continue;}
                        }
                        else if(className == nullptr && !scMustHaveClass)
                        {
                            className = "wildtype";
                        }
                    }

                    //add ABSc dataLine
                    writeToRawDataLock.lock(); //maybe better lock inside the rawData (keep in mind)
                    rawData.add_to_scAbDict(it->second, umiCountMap.at(uniqueID), className);
                    writeToRawDataLock.unlock();

                    readsToKeep += umiCountMap.at(uniqueID);

                    //ONLY the first encountered real read with unique UMI (>90%) is written into dict
                    if(umiFilterThreshold>0.5){break;}
            }
        }
        if(readsWithNoClass > 0)
        {
            result.add_removed_reads_class(readsWithNoClass);
        }
    }
    result.add_removed_reads_umi(totalReadCount - (readsToKeep + readsWithNoClass) );

    if( (totalCount >= 100) && (count % (totalCount / 100) == 0) )
    {
        statusUpdateLock.lock();
        double perc = count/ (double) totalCount;
        printProgress(perc);
        statusUpdateLock.unlock();
    }
    ++count;
}


void BarcodeProcessingHandler::count_umi_occurence(std::vector<int>& positionsOfSameUmi, 
                                                   umiCount& umiLineTmp,
                                                   const std::vector<dataLinePtr>& allScAbCounts,
                                                   const int& lastIdx)
{
    unsigned long long numberAlignedUmis = 0;
    for(int j = 0; j < (allScAbCounts.size() - 1); ++j)
    {
        //calling outputSense algorithm, much faster than levenshtein O(e*max(m,n))
        //however is recently implemented without backtracking
        //before umiMismatches was increased by the length difference between the two UMIs 
        //(no longer done, those deletion should probably be considered as part of the allowed umiMismatches)
        const char* umia = allScAbCounts.at(lastIdx)->umiSeq;
        const char* umib = allScAbCounts.at(j)->umiSeq;
        unsigned int dist = UINT_MAX;
        bool similar = outputSense(umia, umib, barcodeInformation.umiMismatches, dist);

        //if mismatches are within range, change UMI seq
        //the new 'correct' UMI sequence is the one of umiLength, if both r of
        //same length, its the first occuring UMI
        //similar only if dist <= barcodeInformation.umiMismatches
        if(similar)
        {
            if(dist > 0)
            {
                ++numberAlignedUmis;
            }
            //UMIs are not corrected in rawData (the rawData keeps the 'wrong' umi sequences)
            //they could be changed by calling 'changeUmi'
            positionsOfSameUmi.push_back(j);     
            
            umiLineTmp.abCount += allScAbCounts.at(j)->umiCount; // increase count for this UMI
        }
    }
    result.add_umi_mismatches(numberAlignedUmis);
}

void BarcodeProcessingHandler::count_abs_per_single_cell(const std::vector<dataLinePtr>& uniqueAbSc,
                                                        std::atomic<unsigned long long>& count,
                                                        const unsigned long long& totalCount,
                                                        std::shared_ptr<std::unordered_map<const char*, std::vector<umiDataLinePtr>, 
                                                        CharHash, CharPtrComparator>>& umiMap)
{
        //correct for UMI mismatches and fill the AbCountvector
        //iterate through same AbScIdx, calculate levenshtein dist for all UMIs and match those with a certain number of mismatches

        //all dataLines for this AB SC combination
        std::vector<dataLinePtr> scAbCounts = uniqueAbSc;

        //data structures to be filled for the UMI and AB count
        scAbCount abLineTmp; // we fill only this one AB SC count
        umiCount umiLineTmp;

        abLineTmp.scID = umiLineTmp.scID = uniqueAbSc.at(0)->scID;
        abLineTmp.abName = umiLineTmp.abName = uniqueAbSc.at(0)->abName;
        abLineTmp.treatment = umiLineTmp.treatment = uniqueAbSc.at(0)->treatmentName;
        abLineTmp.className = uniqueAbSc.at(0)->cellClassname;

        //if we have no umis erase whole vector and count every element
        if(std::string(scAbCounts.back()->umiSeq) == "" )
        {
            abLineTmp.abCount = scAbCounts.size();
            scAbCounts.clear();
        }
        else
        {
            //sort vector by distance to origional UMI length: 
            //reads get a value for dsitance to UMi length (one Base plus minus gets same value)
            //and are then sorted in decreasing fashion (when comparing UMIs a UMI of length umilength is chosen first)
            //to minimize erros bcs e.g. three reads are within 2MM but we choose one UMI out the outer end regarding MM
            sort(scAbCounts.rbegin(), scAbCounts.rend(), less_than_umi(barcodeInformation.umiLength, umiMap));
        }
        //we take always last element in vector of read of same AB and SC ID
        //then store all reads where UMIs are within distance, and delete those line, and sum up the AB count by one
        while(!scAbCounts.empty())
        {
            dataLinePtr lastAbSc = scAbCounts.back();
            int lastIdx = scAbCounts.size() - 1;

            //check if we have to delete element anyways, bcs umi is too long
            if( std::abs(int( strlen(lastAbSc->umiSeq) - barcodeInformation.umiLength )) >  barcodeInformation.umiMismatches)
            {
                //scAbCounts.pop_back();
                //continue;
            }

            //if in last element
            if(scAbCounts.size() == 1)
            {
                ++abLineTmp.abCount;
                umiLineTmp.abCount += lastAbSc->umiCount; 
                umiLineTmp.umi = lastAbSc->umiSeq;

                result.add_umi_count(umiLineTmp);
                result.add_umi_stats(umiLineTmp);
                scAbCounts.pop_back();
                break;
            }

            //otherwise conmpare all and mark the ones to delete
            std::vector<int> deletePositions;
            umiLineTmp.abCount += lastAbSc->umiCount; //count the first occurence

            //count all occurences of the last UMI for this AB-SC:
            // store the positions of other UMIs within umiMismatches to delete and not count
            //those to final sc AB count
            if(barcodeInformation.umiMismatches > 0)
            {
                count_umi_occurence(deletePositions, umiLineTmp, scAbCounts, lastIdx);
            }
            deletePositions.push_back(lastIdx);

            //ADD UMI if exists
            if(umiLineTmp.abCount > 0)
            {
                umiLineTmp.umi = lastAbSc->umiSeq;
                result.add_umi_count(umiLineTmp);
                result.add_umi_stats(umiLineTmp);
                umiLineTmp.abCount = 0; //reset the count for this umi, all other variables stay the same (Ab name, cellID, etc...)
            }

            //delte all same UMIs
            for(int posIdx = (deletePositions.size() - 1); posIdx >= 0; --posIdx)
            {
                int pos = deletePositions.at(posIdx);
                scAbCounts.erase(scAbCounts.begin() + pos);
            }

            //increase AB count for this one UMI
            ++abLineTmp.abCount;
        }

        //add the data to AB counts if it exists
        if(abLineTmp.abCount>0)
        {
            result.add_ab_count(abLineTmp);
        }

        if((totalCount >= 100) && (count % (totalCount / 100) == 0))
        {
            statusUpdateLock.lock();
            double perc = count/ (double) totalCount;
            printProgress(perc);
            statusUpdateLock.unlock();
        }
        ++count;
}

void BarcodeProcessingHandler::processBarcodeMapping(const int& thread)
{

    //implement thread safe update functions for the data
    //Abcounts Umicounts UmiLog

    //check all UMIs and keep their reads only if they are for >90% a unique scID, ABname, treatmentname
    //also collapse UMIs, and assign the className to each new dataLine, dataLines are now stored as const lines and can no longer be changed
    std::atomic<unsigned long long> umiCount = 0; //using atomic<int> as thread safe read count
    unsigned long long totalCount = rawData.getUniqueUmis()->size();

    boost::asio::thread_pool pool_1(thread); //create thread pool
    std::cout << "STEP[2/3]\t(Remove all reads for a UMI with <90% coming from same AB/SC combination)\n";
    const std::shared_ptr< std::unordered_map<const char*, std::vector<umiDataLinePtr>, CharHash, CharPtrComparator>> umiMap = rawData.getUniqueUmis();
    if(!umiMap->empty())
    {
        for(std::unordered_map<const char*, std::vector<umiDataLinePtr>, CharHash, CharPtrComparator>::const_iterator it = umiMap->begin(); 
        it != umiMap->end(); 
        it++)    
        {
            //careful: uniqueUmi goe sout of scope after enqueuing, therefore just copied...
            boost::asio::post(pool_1, std::bind(&BarcodeProcessingHandler::markReadsWithNoUniqueUmi, this, 
                                            std::cref(it->second), 
                                            std::ref(umiCount), std::cref(totalCount)));
        }
        pool_1.join();
        printProgress(1);
        std::cout << "\n";
    }

    //generate ABcounts per single cell:
    umiCount = 0; //using atomic<int> as thread safe read count
    totalCount = rawData.getUniqueAbSc()->size();
    boost::asio::thread_pool pool_3(thread); //create thread pool
    std::cout << "STEP[3/3]\t(Count reads for AB in single cells)\n";
            
    std::shared_ptr<std::unordered_map<const char*, std::vector<umiDataLinePtr>, 
                    CharHash, CharPtrComparator>> uniqueUmiMap = 
                    rawData.getUniqueUmis();

    const std::shared_ptr< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator>> AbScMap = rawData.getUniqueAbSc();
    for(std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator>::const_iterator it = AbScMap->begin(); 
        it != AbScMap->end(); 
        it++)
    {
        //as above: abSc is copied only
        boost::asio::post(pool_3, std::bind(&BarcodeProcessingHandler::count_abs_per_single_cell, this, 
                                            std::cref(it->second), std::ref(umiCount), std::cref(totalCount), 
                                            std::ref(uniqueUmiMap) ));
    }
    pool_3.join();
    printProgress(1);
    std::cout << "\n";

}

bool BarcodeProcessingHandler::checkIfLineIsDeleted(const dataLinePtr& line, const std::vector<dataLinePtr>& dataLinesToDelete)
{
    if(std::find(dataLinesToDelete.begin(), dataLinesToDelete.end(),line) != dataLinesToDelete.end())
    {
        return true;
    }
    return false;
}

void BarcodeProcessingHandler::writeLog(std::string output)
{
    //WRITE INTO FILE
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");
    if(found == std::string::npos)
    {
        output = "LOG" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "LOG" + output.substr(found+1);
    }
    outputFile.open (output);

    outputFile << "TOTAL READS:\t" << result.get_log_data().totalReads << "\n";
    outputFile << "TOTAL GUIDE-READS:\t" << result.get_log_data().totalGuideReads << "\n";
    outputFile << "TOTAL AB-READS:\t" << result.get_log_data().totalAbReads << "\n";

    outputFile << "REMOVED READS because UMI was not unique(>=90% reads are from same AB/SC):\t" << result.get_log_data().removedUmiReads << "\n";
    outputFile << "REMOVED READS because the single cell had no mapped CLASS:\t" << result.get_log_data().removedClassReads << "\n";

    outputFile << "UMI MM:\t" << result.get_log_data().umiMM << "\n";
    outputFile << "Lost SINGLECELL->CLASS mappings because guide reads for single cell were not unique(>=90% reads for same CLASS):\t" << result.get_log_data().removedClasses << "\n";

    outputFile.close();
}

void BarcodeProcessingHandler::writeAbCountsPerSc(const std::string& output)
{
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");

    //STORE RAW UMI CORRECTED DATA
    std::string umiOutput = output;
    if(found == std::string::npos)
    {
        umiOutput = "UMI" + output;
    }
    else
    {
        umiOutput = output.substr(0,found) + "/" + "UMI" + output.substr(found+1);
    }
    outputFile.open (umiOutput);
    outputFile << "UMI" << "\t" << "AB" << "\t" << "SingleCell_ID" << "\t" << "TREATMENT" << "\t" << "UMI_COUNT" << "\n"; 
    for(umiCount line : result.get_umi_data())
    {
        outputFile << line.umi << "\t" << line.abName << "\t" << line.scID << "\t" << line.treatment << "\t" << line.abCount << "\n"; 
    }
    outputFile.close();

    //STORE AB COUNT DATA
    std::string abOutput = output;
    if(found == std::string::npos)
    {
        abOutput = "AB" + output;
    }
    else
    {
        abOutput = output.substr(0,found) + "/" + "AB" + output.substr(found+1);
    }
    outputFile.open (abOutput);
    bool writeClassLabels = rawData.check_class();
    if(writeClassLabels)
    {
        outputFile << "AB_BARCODE" << "\t" << "SingleCell_BARCODE" << "\t" << "AB_COUNT" << "\t" << "TREATMENT" << "\t" << "CLASS" << "\t" << "CLASS_COUNT" <<"\n"; 
    }
    else
    {
        outputFile << "AB_BARCODE" << "\t" << "SingleCell_BARCODE" << "\t" << "AB_COUNT" << "\t" << "TREATMENT" << "\n"; 
    }
    for(scAbCount line : result.get_ab_data())
    {
        if(writeClassLabels)
        {
                outputFile << line.abName << "\t" << line.scID << "\t" << line.abCount << "\t" << line.treatment << "\t" << line.className << "\t" << guideCountPerSC.at(line.scID) << "\n"; 
        }
        else
        {
                outputFile << line.abName << "\t" << line.scID << "\t" << line.abCount << "\t" << line.treatment << "\n"; 
        }
    }
    outputFile.close();
    
    //store statistics like UMI counts
    std::string umiStatOutput = output;
    if(found == std::string::npos)
    {
        umiOutput = "UMISTAT" + output;
    }
    else
    {
        umiOutput = output.substr(0,found) + "/" + "UMISTAT" + output.substr(found+1);
    }
    outputFile.open (umiOutput);
    outputFile << "UMI_AMPLIFICATION" << "\t" << "AB" << "\t" << "OCCURENCE" << "\n"; 
    umiDist stats = result.get_umi_stats();
    for (auto it : (stats.abs))
    {
        for (auto it2 : stats.dist.at(it.second))
        {
            outputFile << it2.first << "\t" << it.first << "\t" << it2.second << "\n";
        }
    }
    outputFile.close();

}
