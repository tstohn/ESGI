#include "BarcodeProcessingHandler.hpp"

#include <unordered_set>
#include <unordered_map>
#include <set>

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

void generateBarcodeDicts(std::string barcodeFile, std::string barcodeIndices, NBarcodeInformation& barcodeIdData, 
                          std::vector<std::string>& proteinDict, const int& protIdx, 
                          std::vector<std::string>* treatmentDict, const int& treatmentIdx)
{
    //parse barcode file into a vector of a vector of all sequences
    std::vector<std::vector<std::string> > barcodeList;
    std::ifstream barcodeFileStream(barcodeFile);
    for(std::string line; std::getline(barcodeFileStream, line);)
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
        barcodeList.push_back(seqVector);
        seqVector.clear();
    }
    barcodeFileStream.close();

    //parse the indices of CI-barcodes (the index of the lines that store CI-barcodes)
    std::stringstream ss;
    ss.str(barcodeIndices);
    while(ss.good())
    {
        std::string substr;
        getline(ss, substr, ',' );
        barcodeIdData.NBarcodeIndices.push_back(stoi(substr));
    }

    //write a mapping of barcode sequence to a unique number for each
    //CI barcoding round
    for(const int& i : barcodeIdData.NBarcodeIndices)
    {
        int barcodeCount = 0;
        std::unordered_map<std::string, int> barcodeMap;
        //for all options of this barcode
        for(const std::string& barcodeEntry : barcodeList.at(i))
        {
            barcodeMap.insert(std::pair<std::string, int>(barcodeEntry,barcodeCount));
            ++barcodeCount;
        }
        barcodeIdData.barcodeIdDict.push_back(barcodeMap);
    }

    barcodeIdData.NTreatmentIdx = treatmentIdx;
    barcodeIdData.NAbIdx = protIdx;

    proteinDict = barcodeList.at(protIdx);

    if(treatmentIdx!=INT_MAX)
    {
        *treatmentDict = barcodeList.at(treatmentIdx);
    }
}

void BarcodeProcessingHandler::generate_unique_sc_to_class_dict(const std::unordered_map< const char*, 
                                                                std::unordered_map< const char*, unsigned long long>>& scClasseCountDict)
{
    std::unordered_map< const char*, const char*> scClassDict;

    //calculate real class for each single cell
    std::unordered_map< const char*, std::unordered_map< const char*, unsigned long long>>::const_iterator it;
    for (it = scClasseCountDict.begin(); it != scClasseCountDict.end(); it++)
    {
        //count the read for each class for the single cell
        std::unordered_map< const char*, unsigned long long>::const_iterator it2;
        unsigned long long totalReads = 0;
        for(it2 = it->second.begin(); it2 != it->second.end(); it2++)
        {
            totalReads += it2->second;
        }
        //check if any class represents more than 90% of reads
        const char* realClass = nullptr;
        bool foundUniqueClass = false;
        for(it2 = it->second.begin(); it2 != it->second.end(); it2++)
        {
            if( (it2->second/totalReads) >= 0.9 )
            {
                realClass = it2->first;
                foundUniqueClass = true;
            }
        }

        //add class for single cell to dict
        if(foundUniqueClass)
        {
            scClassDict.insert(std::pair<const char*, const char*>(it->first, realClass));
        }
        else
        {
            result.add_removed_class_for_single_cell();
        }
    }

    rawData.set_cell_to_class_dict(scClassDict);

}

void BarcodeProcessingHandler::parseFile(const std::string fileName, const int& thread)
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
    
    std::unordered_map< const char*, std::unordered_map< const char*, unsigned long long>> scClasseCountDict;
    parseBarcodeLines(&instream, totalReads, currentReads, scClasseCountDict);

    //finally add the class of each single cell
    if(rawData.check_class())
    {
        generate_unique_sc_to_class_dict(scClasseCountDict);
    }

    file.close();
}

void BarcodeProcessingHandler::parseBarcodeLines(std::istream* instream, const unsigned long long& totalReads, unsigned long long& currentReads, 
                                                    std::unordered_map< const char*, std::unordered_map< const char*, unsigned long long>>& scClasseCountDict)
{
    std::string line;
    std::cout << "STEP[1/3]\t(READING ALL LINES INTO MEMORY)\n";
    int count = 0;
    int elements = 0; //check that each row has the correct number of barcodes
    unsigned long long abReadCount = 0;
    unsigned long long guideReadCount = 0;
    while(std::getline(*instream, line))
    {
        //for the first read check the positions in the string that refer to CIBarcoding positions
        if(currentReads==0){
            ++currentReads; 
            getBarcodePositions(line, elements);
            continue;
        }
        add_line_to_temporary_data(line, elements, scClasseCountDict, abReadCount, guideReadCount);   

        double perc = currentReads/ (double)totalReads;
        ++currentReads;
        printProgress(perc);        
    }

    result.set_total_reads(currentReads-1); //minus header line
    result.set_total_ab_reads(abReadCount);
    result.set_total_guide_reads(guideReadCount);

    printProgress(1);
    std::cout << "\n";
}

void BarcodeProcessingHandler::add_line_to_temporary_data(const std::string& line, const int& elements,
   std::unordered_map< const char*, std::unordered_map< const char*, unsigned long long>>& scClasseCountDict,
   unsigned long long& abReadCount, unsigned long long& guideReadCount)
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
        std::cerr << "Error in barcode file, following row has not the correct number of sequences: " << line << "\n";
        exit(EXIT_FAILURE);
    }

    //hand over the UMI string, ab string, singleCellstring (concatenation of CIbarcodes)
    std::vector<std::string> ciBarcodes;
    for(int i : fastqReadBarcodeIdx)
    {
        ciBarcodes.push_back(result.at(i));
    }
    std::string singleCellIdx = generateSingleCellIndexFromBarcodes(ciBarcodes);
    
    std::string proteinName = "";
    if(rawData.check_class())
    {
        bool classLine = false;
        std::string name = rawData.get_protein_or_class_name(result.at(abIdx), classLine);
        if(classLine)
        {
            ++guideReadCount;
            rawData.add_tmp_class_line(name, singleCellIdx, scClasseCountDict);
            return;
        }
        else
        {
            proteinName = name;
        }
    }
    else
    {
        proteinName = rawData.getProteinName(result.at(abIdx));
    }
    
    std::string treatment = "";
    if(treatmentIdx != INT_MAX)
    {
        treatment = rawData.getTreatmentName(result.at(treatmentIdx));
    }

    ++abReadCount;
    rawData.add_to_umiDict(result.at(umiIdx), proteinName, singleCellIdx, treatment);
}

std::string BarcodeProcessingHandler::generateSingleCellIndexFromBarcodes(std::vector<std::string> ciBarcodes)
{
    std::string scIdx;

    for(int i = 0; i < ciBarcodes.size(); ++i)
    {
        std::string barcodeAlternative = ciBarcodes.at(i);
        int tmpIdx = varyingBarcodesPos.barcodeIdDict.at(i)[barcodeAlternative];
        scIdx += std::to_string(tmpIdx);
        if(i < ciBarcodes.size() - 1)
        {
            scIdx += ".";
        }
    }

    return scIdx;
}

void BarcodeProcessingHandler::getBarcodePositions(const std::string& line, int& barcodeElements)
{
    std::vector<std::string> result;
    std::stringstream ss;
    ss.str(line);
    int count = 0;
    int variableBarcodeCount = 0;
    std::string substr;
    while(std::getline(ss, substr, '\t'))
    {
        if(substr.empty()){continue;}
        //if substr is only N's
        if(substr.find_first_not_of('N') == std::string::npos)
        {
            //add index for treatment
            if(variableBarcodeCount == varyingBarcodesPos.NTreatmentIdx)
            {
                treatmentIdx = count;
            }
            //add index for combinatorial indexing barcode
            if (std::count(varyingBarcodesPos.NBarcodeIndices.begin(), varyingBarcodesPos.NBarcodeIndices.end(), variableBarcodeCount)) 
            {
                fastqReadBarcodeIdx.push_back(count);
            }
            if(variableBarcodeCount == varyingBarcodesPos.NAbIdx)
            {
                abIdx = count;
            }
            ++variableBarcodeCount;
        }
        else if(substr.find_first_not_of('X') == std::string::npos)
        {
            umiIdx = count;
            umiLength = strlen(substr.c_str());
        }
        result.push_back( substr );
        ++count;
    }
    barcodeElements = count;
    assert(fastqReadBarcodeIdx.size() == varyingBarcodesPos.NBarcodeIndices.size());
    assert(umiIdx != INT_MAX);
    assert(abIdx != INT_MAX);
}

//all reads of the same UMI are combined -> written to a dict for an AB of a unique cell (this read is stored only once, there can already be seen
//as a UMI collapsing step)
void BarcodeProcessingHandler::markReadsWithNoUniqueUmi(const std::vector<umiDataLinePtr>& uniqueUmis,
                                                        std::atomic<unsigned long long>& count,
                                                        const unsigned long long& totalCount)
{
    //count how often we see which AB-SC combinations for this certain UMI
    std::unordered_map<std::string, unsigned long long> umiCountMap;
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
        }
    }

    //check there is a single cell + AB combination representing more than 90% of the UMI reads
    std::string realSingleCellID;
    bool realSingleCellExists = false;
    unsigned long long readsToKeep = 0;
    for(auto singleCellCountPair : umiCountMap)
    {
        double singleCellPerc = (double)singleCellCountPair.second/totalReadCount;
        if( singleCellPerc >= 0.9 )
        {
            realSingleCellID = singleCellCountPair.first;
            realSingleCellExists = true;
            readsToKeep = singleCellCountPair.second;
            break;
        }
    }

    //Collapse UMIs, remove false reads, map a single cell class name to single cells
    unsigned long long readsWithNoClass = 0;
    if(realSingleCellExists)
    {
        //delete only the <=10% 'false' reads
        for(int i = 0; i < uniqueUmis.size(); ++i)
        {
            std::string uniqueID = std::string(uniqueUmis.at(i)->scID) + std::string(uniqueUmis.at(i)->abName);
            if(uniqueID == realSingleCellID)      //add to new DIct
            {

                    const char* className = nullptr;
                    if(rawData.check_class())
                    {
                        className = rawData.get_sc_class_name(uniqueUmis.at(i)->scID);
                        if(className == nullptr)
                        {
                            //single cell has no class name
                            readsWithNoClass = readsToKeep;
                            break;
                        }
                    }

                    //add ABSc dataLine
                    writeToRawDataLock.lock(); //maybe better lock inside the rawData (keep in mind)
                    rawData.add_to_scAbDict(uniqueUmis.at(i), readsToKeep, className);
                    writeToRawDataLock.unlock();

                    //ONLY the first encountered real read with unique UMI (>90%) is written into dict
                    break;
            }
        }
        if(readsWithNoClass > 0)
        {
            result.add_removed_reads_class(readsWithNoClass);
        }
    }

    result.add_removed_reads_umi(totalReadCount - readsToKeep);

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
                                                   const int& umiMismatches,
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
        int dist = INT_MAX;
        int start = 0;
        int end = 0;
        bool similar = outputSense(umia, umib, umiMismatches, dist);

        //if mismatches are within range, change UMI seq
        //the new 'correct' UMI sequence is the one of umiLength, if both r of
        //same length, its the first occuring UMI
        if(dist <= umiMismatches)
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

void BarcodeProcessingHandler::count_abs_per_single_cell(const int& umiMismatches, const std::vector<dataLinePtr>& uniqueAbSc,
                                                        std::atomic<unsigned long long>& count,
                                                        const unsigned long long& totalCount,
                                                        std::shared_ptr<std::unordered_map<const char*, std::vector<umiDataLinePtr>, 
                                                        CharHash, CharPtrComparator>>& umiMap)
{
   //correct for UMI mismatches and fill the AbCountvector
    //iterate through same AbScIdx, calculate levenshtein dist for all UMIs and match those with a certain number of mismatches

        //all dataLines for this AB SC combination
        std::vector<dataLinePtr> scAbCounts = uniqueAbSc;

        //sort vector by distance to origional UMI length: 
        //reads get a value for dsitance to UMi length (one Base plus minus gets same value)
        //and are then sorted in decreasing fashion (when comparing UMIs a UMI of length umilength is chosen first)
        //to minimize erros bcs e.g. three reads are within 2MM but we choose one UMI out the outer end regarding MM
        sort(scAbCounts.rbegin(), scAbCounts.rend(), less_than_umi(umiLength, umiMap));

        //data structures to be filled for the UMI and AB count
        scAbCount abLineTmp; // we fill only this one AB SC count
        umiCount umiLineTmp;

        abLineTmp.scID = umiLineTmp.scID = uniqueAbSc.at(0)->scID;
        abLineTmp.abName = umiLineTmp.abName = uniqueAbSc.at(0)->abName;
        abLineTmp.treatment = umiLineTmp.treatment = uniqueAbSc.at(0)->treatmentName;
        abLineTmp.className = uniqueAbSc.at(0)->cellClassname;

        //we take always last element in vector of read of same AB and SC ID
        //then store all reads wwhere UMIs are within distance, and delete those line, and sum up the AB count by one
        while(!scAbCounts.empty())
        {

            dataLinePtr lastAbSc = scAbCounts.back();
            int lastIdx = scAbCounts.size() - 1;
            //check if we have to delete element anyways, since it is not a unique UMI
            //in this case we simply delete this line and check the next one

            //if in last element
            if(scAbCounts.size() == 1)
            {
                ++abLineTmp.abCount;
                umiLineTmp.abCount += lastAbSc->umiCount; 
                umiLineTmp.umi = lastAbSc->umiSeq;

                result.add_umi_count(umiLineTmp);
                scAbCounts.pop_back();
                break;
            }

            //otherwise conmpare all and mark the ones to delete
            std::vector<int> deletePositions;
            umiLineTmp.abCount += lastAbSc->umiCount; //count the first occurence
            //count all occurences of the last UMI for this AB-SC
            count_umi_occurence(deletePositions, umiLineTmp, scAbCounts, umiMismatches, lastIdx);
            deletePositions.push_back(lastIdx);

            //ADD UMI if exists
            if(umiLineTmp.abCount > 0)
            {
                umiLineTmp.umi = lastAbSc->umiSeq;
                result.add_umi_count(umiLineTmp);
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

void BarcodeProcessingHandler::processBarcodeMapping(const int& umiMismatches, const int& thread)
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
    for(std::unordered_map<const char*, std::vector<umiDataLinePtr>, CharHash, CharPtrComparator>::const_iterator it = umiMap->begin(); 
        it != umiMap->end(); 
        it++)    {
        //careful: uniqueUmi goe sout of scope after enqueuing, therefore just copied...
        boost::asio::post(pool_1, std::bind(&BarcodeProcessingHandler::markReadsWithNoUniqueUmi, this, 
                                          std::cref(it->second), 
                                          std::ref(umiCount), std::cref(totalCount)));
    }
    pool_1.join();
    printProgress(1);
    std::cout << "\n";

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
                                          std::ref(umiMismatches), std::cref(it->second), 
                                          std::ref(umiCount), std::cref(totalCount), std::ref(uniqueUmiMap) ));
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
        outputFile << "AB_BARCODE" << "\t" << "SingleCell_BARCODE" << "\t" << "AB_COUNT" << "\t" << "TREATMENT" << "\t" << "CLASS" << "\n"; 
    }
    else
    {
        outputFile << "AB_BARCODE" << "\t" << "SingleCell_BARCODE" << "\t" << "AB_COUNT" << "\t" << "TREATMENT" << "\n"; 
    }

    for(scAbCount line : result.get_ab_data())
    {
        if(writeClassLabels)
        {
                outputFile << line.abName << "\t" << line.scID << "\t" << line.abCount << "\t" << line.treatment << "\t" << line.className << "\n"; 
        }
        else
        {
                outputFile << line.abName << "\t" << line.scID << "\t" << line.abCount << "\t" << line.treatment << "\n"; 
        }
    }
    outputFile.close();
}