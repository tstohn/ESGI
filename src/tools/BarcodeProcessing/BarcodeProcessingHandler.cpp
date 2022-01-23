#include "BarcodeProcessingHandler.hpp"
#include "helper.hpp"

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

void writeExtendedData(const std::string& output, std::vector<std::vector<std::pair<unsigned long long, unsigned long long>>>& uniqueUmiToAbScVec, 
                       std::vector<std::vector<umiQualityExtended>> extendedQualityVec)
{
    //write everything to a file
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");
    //STORE RAW UMI CORRECTED DATA
    //TODO: delete in beginning if already there !!!
    //STORE BARCODE SPECIFIC
    std::string umiOutput1 = output;
    if(found == std::string::npos)
    {
        umiOutput1 = "UMIQUAL_" + output;
    }
    else
    {
        umiOutput1 = output.substr(0,found) + "/" + "UMIQUAL_" + output.substr(found+1);
    }
    std::remove(umiOutput1.c_str());
    outputFile.open (umiOutput1, std::ofstream::app);
    outputFile << "NumberOfReadsForUMI\tNumberOfDifferentAbSc\n";
    for(int i = 0; i < uniqueUmiToAbScVec.size(); ++i)
    {
        for(int j = 0; j < uniqueUmiToAbScVec.at(i).size(); ++j)
        {
            outputFile << uniqueUmiToAbScVec.at(i).at(j).first << "\t" << uniqueUmiToAbScVec.at(i).at(j).second << "\n"; 
        }
    }
    outputFile.close();

    //STORE UMI COUNTS
    std::string umiOutput2 = output;
    if(found == std::string::npos)
    {
        umiOutput2 = "UMICOUNT_" + output;
    }
    else
    {
        umiOutput2 = output.substr(0,found) + "/" + "UMICOUNT_" + output.substr(found+1);
    }
    std::remove(umiOutput2.c_str());
    outputFile.open (umiOutput2, std::ofstream::app);
    outputFile << "UmiNumber\tBCAb\tBC4\tBC1\tBC2\tBC3\tBCAb_30\tBCAb_60\tBCAb_90\tBC4_30\tBC4_60\tBC4_90\tBC1_30\tBC1_60\tBC1_90\tBC2_30\tBC2_60\tBC2_90\tBC3_30\tBC3_60\tBC3_90\n";
    for(int i = 0; i < uniqueUmiToAbScVec.size(); ++i)
    {
        for(int j = 0; j < uniqueUmiToAbScVec.at(i).size(); ++j)
        {
            umiQualityExtended qualEx = extendedQualityVec.at(i).at(j);
            outputFile << qualEx.numOfUmiOccurences << "\t" << qualEx.numOfABDifferences << "\t" << 
            qualEx.numOfBC4Differences << "\t" << qualEx.numOfBC1Differences << "\t" << 
            qualEx.numOfBC2Differences << "\t" << qualEx.numOfBC3Differences << "\t";

            // add a list with how many different barcodes you have to add to reach 30%, 60%, 90% of data
            //meaning: how many percent of the number of different barcodes do i need, to get 30,etc. % of the number of reads 
            outputFile << calcualtePercentages(qualEx.AbBarcodeDistribution, qualEx.numOfUmiOccurences, 0.3) << "\t";
            outputFile << calcualtePercentages(qualEx.AbBarcodeDistribution, qualEx.numOfUmiOccurences, 0.6) << "\t";
            outputFile << calcualtePercentages(qualEx.AbBarcodeDistribution, qualEx.numOfUmiOccurences, 0.9) << "\t";

            outputFile << calcualtePercentages(qualEx.BC4Distribution, qualEx.numOfUmiOccurences, 0.3) << "\t";
            outputFile << calcualtePercentages(qualEx.BC4Distribution, qualEx.numOfUmiOccurences, 0.6) << "\t";
            outputFile << calcualtePercentages(qualEx.BC4Distribution, qualEx.numOfUmiOccurences, 0.9) << "\t";

            outputFile << calcualtePercentages(qualEx.BC1Distribution, qualEx.numOfUmiOccurences, 0.3) << "\t";
            outputFile << calcualtePercentages(qualEx.BC1Distribution, qualEx.numOfUmiOccurences, 0.6) << "\t";
            outputFile << calcualtePercentages(qualEx.BC1Distribution, qualEx.numOfUmiOccurences, 0.9) << "\t";

            outputFile << calcualtePercentages(qualEx.BC2Distribution, qualEx.numOfUmiOccurences, 0.3) << "\t";
            outputFile << calcualtePercentages(qualEx.BC2Distribution, qualEx.numOfUmiOccurences, 0.6) << "\t";
            outputFile << calcualtePercentages(qualEx.BC2Distribution, qualEx.numOfUmiOccurences, 0.9) << "\t";

            outputFile << calcualtePercentages(qualEx.BC3Distribution, qualEx.numOfUmiOccurences, 0.3) << "\t";
            outputFile << calcualtePercentages(qualEx.BC3Distribution, qualEx.numOfUmiOccurences, 0.6) << "\t";
            outputFile << calcualtePercentages(qualEx.BC3Distribution, qualEx.numOfUmiOccurences, 0.9) << "\n";
        }
    }
    outputFile.close();
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

    if(treatmentDict != nullptr)
    {
        *treatmentDict = barcodeList.at(treatmentIdx);
    }

}

void BarcodeProcessingHandler::parseFile(const std::string fileName, const int& thread)
{
    int totalReads = totalNumberOfLines(fileName);
    int currentReads = 0;
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

void BarcodeProcessingHandler::parseBarcodeLines(std::istream* instream, const int& totalReads, int& currentReads)
{
    std::string line;
    std::cout << "READING ALL LINES INTO MEMORY\n";
    int count = 0;
    int elements = 0; //check that each row has the correct number of barcodes
    while(std::getline(*instream, line))
    {
        //for the first read check the positions in the string that refer to CIBarcoding positions
        if(currentReads==0){
            ++currentReads; 
            getBarcodePositions(line, elements);
            continue;
        }
        addFastqReadToUmiData(line, elements);   

        double perc = currentReads/ (double)totalReads;
        ++currentReads;
        printProgress(perc);        
    }
    printProgress(1);
    std::cout << "\n";
}

void BarcodeProcessingHandler::addFastqReadToUmiData(const std::string& line, const int& elements)
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
    std::string proteinName = rawData.getProteinName(result.at(abIdx));
    
    std::string treatment = "";

    if(treatmentIdx != INT_MAX)
    {
        treatment = rawData.getTreatmentName(result.at(treatmentIdx));
    }
    rawData.add(result.at(umiIdx), proteinName, singleCellIdx, treatment);

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

void BarcodeProcessingHandler::processBarcodeMapping(const int& umiMismatches, const int& thread)
{
    //ToDO: can be done with a thread queue
    //split AbSc map into equal parts
    std::vector< std::vector < std::vector<dataLinePtr> > > independantAbScBatches(thread);
    std::vector< std::vector < std::vector<dataLinePtr> > > independantUmiBatches(thread);

    int element_number = rawData.getUniqueAbSc().size() / thread;

    //iterate through map, and add one element after the other ot the vector (each element is itself a vector
    //of all the lines with the same AB and SingleCell)    
    int abScLineIdx = 0;
    for(auto mapElement : rawData.getUniqueAbSc())
    {
        std::vector<dataLinePtr> abScLine = mapElement.second;
        
        independantAbScBatches.at( abScLineIdx%thread ).push_back(abScLine);

        ++abScLineIdx;
    }

    int umiLineIdx = 0;
    for(auto mapElement : rawData.getUniqueUmis())
    {
        std::vector<dataLinePtr> umiLine = mapElement.second;
        
        independantUmiBatches.at( umiLineIdx%thread ).push_back(umiLine);

        ++umiLineIdx;
    }

    //temporary vectors to store all the thread outputs
    std::vector<StatsUmi> umiStatsThreaded(thread);
    std::vector<umiQuality> umiQualThreaded(thread);

    std::vector<std::vector<dataLinePtr> > umiDataThreaded(thread);
    std::vector<std::vector<scAbCount> > abDataThreaded(thread);

    std::vector<std::vector<dataLinePtr> > dataLinesToDeleteThreaded(thread);

    //for every batch calcualte the unique umis with same/different AbSc barcodes
    std::vector<std::thread> workers;
    int currentUmisChecked = 0;
    std::cout << "Checking Quality of UMIs\n";
    for (int i = 0; i < thread; ++i) 
    {
        workers.push_back(std::thread(&BarcodeProcessingHandler::umiQualityCheck, this, std::ref(independantUmiBatches.at(i)), std::ref(umiQualThreaded.at(i)), std::ref(currentUmisChecked) ));
    }
    printProgress(1);
    for (std::thread &t: workers) 
    {
        if (t.joinable()) {
            t.join();
        }
    }
    std::cout << "\n";

    //get all datalinePtrs that should be deleted since they are probably background noise (additional reads that are not 90% of AbSC combination for a UMI)
    std::cout << "Removing false reads by wrong Single Cells IDs for reads of same UMI\n";
    int currentReadsRemoved = 0;
    std::vector<dataLinePtr> dataLinesToDelete;
    for (int i = 0; i < thread; ++i) 
    {
        workers.push_back(std::thread(&BarcodeProcessingHandler::removeFalseSingleCellsFromUmis, this, std::ref(independantUmiBatches.at(i)), std::ref(currentReadsRemoved),
        std::ref(dataLinesToDeleteThreaded.at(i))
        ));
    }
    printProgress(1);
    for (std::thread &t: workers) 
    {
        if (t.joinable()) {
            t.join();
        }
    }
    for (int i = 0; i < thread; ++i) 
    {
        //combine threaded data
        dataLinesToDelete.insert(dataLinesToDelete.end(), dataLinesToDeleteThreaded.at(i).begin(), dataLinesToDeleteThreaded.at(i).end());
    }
    std::cout << "\n";


    //for every batch calculate an UnprocessedDemultiplexedData and ABData vector and stats
    int currentUmisCorrected = 0;
    std::cout << "Correcting UMIs and counting ABs\n";
    for (int i = 0; i < thread; ++i) 
    {
        workers.push_back(std::thread(&BarcodeProcessingHandler::correctUmis, this, std::ref(umiMismatches), std::ref(umiStatsThreaded.at(i)), std::ref(umiDataThreaded.at(i)),
                          std::ref(abDataThreaded.at(i)), std::ref(independantAbScBatches.at(i)), std::ref(currentUmisCorrected), std::ref(dataLinesToDelete)));
    }
    printProgress(1);
    for (std::thread &t: workers) 
    {
        if (t.joinable()) {
            t.join();
        }
    }
    std::cout << "\n";

    //combine the three dataSets
    for (int i = 0; i < thread; ++i) 
    {
        //data
        umiData.insert(umiData.end(), umiDataThreaded.at(i).begin(), umiDataThreaded.at(i).end());
        abData.insert(abData.end(), abDataThreaded.at(i).begin(), abDataThreaded.at(i).end());

        //statistics
        //simply add integer values
        ullong_save_add(qual.sameUmiDiffAbSc, umiQualThreaded.at(i).sameUmiDiffAbSc);
        ullong_save_add(qual.sameUmiSameAbSc, umiQualThreaded.at(i).sameUmiSameAbSc);
        //add each element of Dict to master dict
        for(auto dictElem : umiStatsThreaded.at(i).umiMismatchDict)
        {
            if(stats.umiMismatchDict.find(dictElem.first) != stats.umiMismatchDict.end())
            {
                ullong_save_add(stats.umiMismatchDict[dictElem.first], dictElem.second);
            }
            else
            {
                stats.umiMismatchDict.insert(dictElem);
            }
        }
    }
}

bool BarcodeProcessingHandler::checkDataLineValidityDueToUmiBackground(const dataLinePtr& line, const std::vector<dataLinePtr>& dataLinesToDelete)
{
    if(std::find(dataLinesToDelete.begin(), dataLinesToDelete.end(),line) != dataLinesToDelete.end())
    {
        return true;
    }
    return false;
}

//correct for mismatches in the UMI
void BarcodeProcessingHandler::correctUmis(const int& umiMismatches, StatsUmi& statsTmp, std::vector<dataLinePtr>& umiDataTmp, std::vector<scAbCount>& abDataTmp, 
                                const std::vector<std::vector<dataLinePtr> >& AbScBucket, int& currentUmisCorrected, 
                                const std::vector<dataLinePtr>& dataLinesToDelete)
{
    //correct for UMI mismatches and fill the AbCountvector
    //iterate through same AbScIdx, calculate levenshtein dist for all UMIs and match those with a certain number of mismatches
    int tmpCurrentUmisCorrected = 0;
    int containerSize = rawData.getUniqueAbSc().size();
    for (auto uniqueAbSc : AbScBucket)
    {
        scAbCount abLineTmp;
        abLineTmp.scID = uniqueAbSc.at(0)->scID;
        abLineTmp.abName = uniqueAbSc.at(0)->abName;
        abLineTmp.treatment = uniqueAbSc.at(0)->treatmentName;

        int abCount = 0; // abCount is calculated for every umi one after the other, if an umi is unique, the count is incremented
        //for umis with several occurences the last occurence will increment the count, the very last UMI is never checked and is always
        //incrementing the count, therefore initialized with 1
        
        //we take always first element in vector of read of same AB and SC ID
        //then store all reads wwhere UMIs are within distance, and delete those line, and sum up the AB count by one
        while(!uniqueAbSc.empty())
        {
            //check if we have to delete element anyways, since it s a read of false UMI background
            //in this case we simply delete this line and check the next one
            if(checkDataLineValidityDueToUmiBackground(uniqueAbSc.at(0), dataLinesToDelete))
            {
                uniqueAbSc.pop_back();
                continue;
            }

            //if in last element
            if(uniqueAbSc.size() == 1)
            {
                ++abCount;
                umiDataTmp.push_back(uniqueAbSc.front());
                uniqueAbSc.pop_back();
                break;
            }
            //otherwise conmpare all and mark the ones to delete
            std::vector<int> deletePositions;
            int i =0;
            deletePositions.push_back(i);
            for(int j = 1; j < uniqueAbSc.size(); ++j)
            {
                const char* umia = uniqueAbSc.at(i)->umiSeq;
                const char* umib = uniqueAbSc.at(j)->umiSeq;
                int dist = INT_MAX;
                int start = 0;
                int end = 0;

                const int diff = umiLength - MIN(std::strlen(umia), std::strlen(umib));
                const int lengthCorrectedMismatches = umiMismatches + (sqrt(diff*diff));
                //calling outputSense algorithm, much faster than levenshtein O(e*max(m,n))
                //however is recently implemented without backtracking
                bool similar = outputSense(umia, umib, lengthCorrectedMismatches, dist);

                //if mismatches are within range, change UMI seq
                //the new 'correct' UMI sequence is the one of umiLength, if both r of
                //same length, its the first occuring UMI
                if(dist <= lengthCorrectedMismatches)
                {
                    if(dist!=0)
                    {
                        //get real UMI
                        const char* realUmi;
                        if(strlen(umia)==umiLength)
                        {
                            realUmi = umia;
                            uniqueAbSc.at(j)->umiSeq = uniqueAbSc.at(i)->umiSeq;
                            rawData.changeUmi(umib, umia, uniqueAbSc.at(i));
                        }
                        else if(strlen(umib)==umiLength)
                        {
                            realUmi = umib;
                            uniqueAbSc.at(i)->umiSeq = uniqueAbSc.at(j)->umiSeq;
                            rawData.changeUmi(umia, umib, uniqueAbSc.at(j));
                        }
                        else
                        {
                            realUmi = umia;
                            uniqueAbSc.at(j)->umiSeq = uniqueAbSc.at(i)->umiSeq;
                            rawData.changeUmi(umib, umia, uniqueAbSc.at(i));
                        }   
                    }  
                    deletePositions.push_back(j);               
                }

                //dist is set inside levenshtein, if its <= mismatches=2
                if(statsTmp.umiMismatchDict.find(dist) != statsTmp.umiMismatchDict.end())
                {
                    ullong_save_add(statsTmp.umiMismatchDict.at(dist), 1);
                }
                else
                {
                    statsTmp.umiMismatchDict.insert(std::make_pair(dist, 1));
                }
            }

            for(int posIdx = deletePositions.size()-1; posIdx >=0; --posIdx)
            {
                int pos = deletePositions.at(posIdx);
                umiDataTmp.push_back(uniqueAbSc.at(pos));
                uniqueAbSc.erase(uniqueAbSc.begin()+pos);
            }

            ++abCount;

        }

        abLineTmp.abCount = abCount;
        if(abCount>0){abDataTmp.push_back(abLineTmp);}

        ++tmpCurrentUmisCorrected;
        if((containerSize >= 100) && (tmpCurrentUmisCorrected % (containerSize / 100) == 0))
        {
            lock.lock();
            currentUmisCorrected += tmpCurrentUmisCorrected;
            tmpCurrentUmisCorrected = 0;
            double perc = currentUmisCorrected/ (double) containerSize;
            printProgress(perc);
            lock.unlock();
        }
    }
}

//unused function, origionally used to compare all UMIs with each other additionally to UMI comparison, to generate
//a statistic of UMI distances
/*
void UmiDataParser::correctUmisWithStats(const int& umiMismatches, StatsUmi& statsTmp, std::vector<dataLinePtr>& umiDataTmp, std::vector<scAbCount>& abDataTmp, 
                                const std::vector<std::vector<dataLinePtr> >& AbScBucket, int& currentUmisCorrected)
{
    //correct for UMI mismatches and fill the AbCountvector
    //iterate through same AbScIdx, calculate levenshtein dist for all UMIs and match those with a certain number of mismatches
    int tmpCurrentUmisCorrected = 0;
    int containerSize = rawData.getUniqueAbSc().size();
    for (auto uniqueAbSc : AbScBucket)
    {
        scAbCount abLineTmp;
        abLineTmp.scID = uniqueAbSc.at(0)->scID;
        abLineTmp.abName = rawData.getProteinName(uniqueAbSc.at(0)->abName);
        abLineTmp.treatment = rawData.getTreatmentName(uniqueAbSc.at(0)->treatmentName);

        int abCount = 1; // abCount is calculated for every umi one after the other, if an umi is unique, the count is incremented
        //for umis with several occurences the last occurence will increment the count, the very last UMI is never checked and is always
        //incrementing the count, therefore initialized with 1
        for(int i = 0; i < (uniqueAbSc.size() - 1); ++i)
        {
            //assert(abLineTmp.abName == uniqueAbSc.second.at(i)->abName);
            //assert(abLineTmp.scID == uniqueAbSc.second.at(i)->scID);
            bool unique = true;
            for(int j = i+1; j < uniqueAbSc.size(); ++j)
            {
                const char* umia = uniqueAbSc.at(i)->umiSeq;
                const char* umib = uniqueAbSc.at(j)->umiSeq;

                int dist = INT_MAX;
                int start = 0;
                int end = 0;

                bool similar = levenshtein(umia, umib, umiLength, start, end, dist, true);

                //if(std::strcmp(abLineTmp.scID, "10295") == 0 & std::strcmp((abLineTmp.abName)->c_str(), "CTD1")==0){std::cout << umia << " " << umib << " => "<< dist << " " << start << " " << end <<  "\n";}

                //if mismatches are within range, change UMI seq
                //the new 'correct' UMI sequence is the one of umiLength, if both r of
                //same length, its the first occuring UMI
                if(dist <= umiMismatches)
                {
                    unique = false;
                    if(dist!=0)
                    {
                        //get real UMI
                        const char* realUmi;
                        if(strlen(umia)==umiLength)
                        {
                            realUmi = umia;
                            uniqueAbSc.at(j)->umiSeq = uniqueAbSc.at(i)->umiSeq;
                            rawData.changeUmi(umib, umia, uniqueAbSc.at(i));
                        }
                        else if(strlen(umib)==umiLength)
                        {
                            realUmi = umib;
                            uniqueAbSc.at(i)->umiSeq = uniqueAbSc.at(j)->umiSeq;
                            rawData.changeUmi(umia, umib, uniqueAbSc.at(j));
                        }
                        else
                        {
                            realUmi = umia;
                            uniqueAbSc.at(j)->umiSeq = uniqueAbSc.at(i)->umiSeq;
                            rawData.changeUmi(umib, umia, uniqueAbSc.at(i));
                        }   
                    }                 
                }

                //dist is set inside levenshtein, if its <= mismatches=2
                if(statsTmp.umiMismatchDict.find(dist) != statsTmp.umiMismatchDict.end())
                {
                    ullong_save_add(statsTmp.umiMismatchDict.at(dist), 1);
                }
                else
                {
                    statsTmp.umiMismatchDict.insert(std::make_pair(dist, 1));
                }
            }
            //if(std::strcmp(abLineTmp.scID, "10295") == 0 & std::strcmp((abLineTmp.abName)->c_str(), "CTD1")==0){std::cout << "\n";}
            if(unique)
            {
                ++abCount;
                            //if(std::strcmp(abLineTmp.scID, "10295") == 0 & std::strcmp((abLineTmp.abName)->c_str(), "CTD1")==0){std::cout << "      "<< abCount << "      ADDED\n";;}
            }
            umiDataTmp.push_back(uniqueAbSc.at(i));
        }    
        abLineTmp.abCount = abCount;
        abDataTmp.push_back(abLineTmp);
        umiDataTmp.push_back(uniqueAbSc.at(uniqueAbSc.size() - 1));

        ++tmpCurrentUmisCorrected;
        if( (containerSize >= 100) && (tmpCurrentUmisCorrected % ( containerSize / 100) == 0) )
        {
            lock.lock();
            currentUmisCorrected += tmpCurrentUmisCorrected;
            tmpCurrentUmisCorrected = 0;
            double perc = currentUmisCorrected/ (double) containerSize;
            printProgress(perc);
            lock.unlock();
        }
    }
}*/

void BarcodeProcessingHandler::extended_umi_quality_check(const int& thread, const std::string& output)
{
    //split AbSc map into equal parts
    std::vector< std::vector < std::vector<dataLinePtr> > > independantUmiBatches(thread);

    int element_number = rawData.getUniqueAbSc().size() / thread;

    //iterate through map, and add one element after the other ot the vector (each element is itself a vector
    //of all the lines with the same AB and SingleCell)    
    int umiLineIdx = 0;
    for(auto mapElement : rawData.getUniqueUmis())
    {
        std::vector<dataLinePtr> umiLine = mapElement.second;
        independantUmiBatches.at( umiLineIdx%thread ).push_back(umiLine);
        ++umiLineIdx;
    }

    //vectors to be filled in threads
    std::vector<std::vector<std::pair<unsigned long long, unsigned long long>>> uniqueUmiToDiffAbScVec(thread);
    std::vector<std::vector<umiQualityExtended>> extendedQualityVec(thread);

    //for every batch calcualte the unique umis with same/different AbSc barcodes
    std::vector<std::thread> workers;
    int currentUmisChecked = 0;
    std::cout << "Checking Extended Quality of UMIs\n";
    for (int i = 0; i < thread; ++i) 
    {
        workers.push_back(std::thread(&BarcodeProcessingHandler::umiQualityCheckExtended, this, std::ref(independantUmiBatches.at(i)), 
                          std::ref(currentUmisChecked),
                          std::ref(uniqueUmiToDiffAbScVec.at(i)), std::ref(extendedQualityVec.at(i)), output ));
    }
    printProgress(1);
    for (std::thread &t: workers) 
    {
        if (t.joinable()) {
            t.join();
        }
    }
    std::cout << "\n";

    //combine vector for threads while writing
    writeExtendedData(output, uniqueUmiToDiffAbScVec, extendedQualityVec);

}

void BarcodeProcessingHandler::umiQualityCheckExtended(const std::vector< std::vector<dataLinePtr> >& uniqueUmis,
                                            int& currentUmisChecked, std::vector<std::pair<unsigned long long, unsigned long long>>& uniqueUmiToDiffAbSc, 
                                            std::vector<umiQualityExtended>& extendedQuality, const std::string& output)
{
    //first quality check, does a unique umi have always the same AbScIdx
    int tmpCurrentUmisChecked =0;
    int containerSize = rawData.getUniqueUmis().size();
    for(auto uniqueUmi : uniqueUmis)
    {
        //structures to fill and push back into result vectors
        umiQualityExtended qualExTmp;

        //tmp values during iteration over reads for same UMI
        std::vector<std::string> uniqueUmis_1;
        std::vector<std::string> uniqueUmis_2;
        std::vector<std::string> uniqueUmis_3;
        std::vector<std::string> uniqueUmis_4;
        std::vector<const char*> uniqueUmis_AB;

        std::vector<unsigned long long> AbBarcodeDistribution;
        std::vector<unsigned long long> BC4Distribution;
        std::vector<unsigned long long> BC1Distribution;
        std::vector<unsigned long long> BC2Distribution;
        std::vector<unsigned long long> BC3Distribution;

        std::unordered_set<std::string> AbScVector;
        for(int i = 0; i < uniqueUmi.size(); ++i)
        {

            dataLinePtr a = uniqueUmi.at(i);
            //keep unique absc in a vector
            std::string ab = a->abName;
            std::string absc = ab + a->scID;
            //make a unique set of all AbSc possibilities for all reads of same UMI
            if(AbScVector.find(absc) == AbScVector.end())
            {
                AbScVector.insert(absc);
            }

            //same for all different barcodes
            std::vector<const char*>::iterator it = std::find(uniqueUmis_AB.begin(), uniqueUmis_AB.end(), a->abName);
            if(it == uniqueUmis_AB.end())
            {
                uniqueUmis_AB.push_back(a->abName);
                AbBarcodeDistribution.push_back(1);
            }
            else
            {
                ++AbBarcodeDistribution.at(std::distance(uniqueUmis_AB.begin(), it));
            }

            std::vector<std::string> aVec = splitByDelimiter(a->scID, ".");
            assert(aVec.size() == 4);
            std::vector<std::string>::iterator it2;

            it2 = std::find(uniqueUmis_4.begin(), uniqueUmis_4.end(), aVec.at(0));
            if(it2 == uniqueUmis_4.end())
            {
                uniqueUmis_4.push_back(aVec.at(0));
                BC4Distribution.push_back(1);
            }
            else
            {
                ++BC4Distribution.at(std::distance(uniqueUmis_4.begin(), it2));
            }

            it2 = std::find(uniqueUmis_1.begin(), uniqueUmis_1.end(), aVec.at(1));
            if(it2 == uniqueUmis_1.end())
            {
                uniqueUmis_1.push_back(aVec.at(1));
                BC1Distribution.push_back(1);
            }
            else
            {
                ++BC1Distribution.at(std::distance(uniqueUmis_1.begin(), it2));
            }

            it2 = std::find(uniqueUmis_2.begin(), uniqueUmis_2.end(), aVec.at(2));
            if(it2 == uniqueUmis_2.end())
            {
                uniqueUmis_2.push_back(aVec.at(2));
                BC2Distribution.push_back(1);
            }
            else
            {
                ++BC2Distribution.at(std::distance(uniqueUmis_2.begin(), it2));
            }
                        
            it2 = std::find(uniqueUmis_3.begin(), uniqueUmis_3.end(), aVec.at(3));
            if(it2 == uniqueUmis_3.end())
            {
                uniqueUmis_3.push_back(aVec.at(3));
                BC3Distribution.push_back(1);
            }
            else
            {
                ++BC3Distribution.at(std::distance(uniqueUmis_3.begin(), it2));
            }
        }

        qualExTmp.numOfUmiOccurences = uniqueUmi.size();
        qualExTmp.numOfABDifferences = uniqueUmis_AB.size();
        qualExTmp.numOfBC4Differences = uniqueUmis_4.size();
        qualExTmp.numOfBC1Differences = uniqueUmis_1.size();
        qualExTmp.numOfBC2Differences = uniqueUmis_2.size();
        qualExTmp.numOfBC3Differences = uniqueUmis_3.size();

        qualExTmp.AbBarcodeDistribution = AbBarcodeDistribution;
        qualExTmp.BC4Distribution = BC4Distribution;
        qualExTmp.BC1Distribution = BC1Distribution;
        qualExTmp.BC2Distribution = BC2Distribution;
        qualExTmp.BC3Distribution = BC3Distribution;

        extendedQuality.push_back(qualExTmp);
        uniqueUmiToDiffAbSc.push_back(std::make_pair<unsigned long long, unsigned long long>(uniqueUmi.size(), AbScVector.size()));

        //Update Progress
        ++tmpCurrentUmisChecked;
        if( (containerSize >= 100) && (tmpCurrentUmisChecked % (containerSize / 100) == 0) )
        {
            lock.lock();
            currentUmisChecked += tmpCurrentUmisChecked;
            tmpCurrentUmisChecked = 0;
            double perc = currentUmisChecked/ (double) containerSize;
            printProgress(perc);
            lock.unlock();
        }
    } // UMI end
}

void BarcodeProcessingHandler::umiQualityCheck(const std::vector< std::vector<dataLinePtr> >& uniqueUmis, umiQuality& qualTmp, int& currentUmisChecked)
{
    //first quality check, does a unique umi have always the same AbScIdx
    int tmpCurrentUmisChecked =0;
    int containerSize = rawData.getUniqueUmis().size();
    for(auto uniqueUmi : uniqueUmis)
    {
        //for all AbSc combinations of this unique UMI
        for(int i = 0; i < (uniqueUmi.size() - 1); ++i)
        {
            for(int j = i+1; j < uniqueUmi.size(); ++j)
            {
                dataLinePtr a = uniqueUmi.at(i);
                dataLinePtr b = uniqueUmi.at(j);
                if(strcmp(a->abName, b->abName)==0 && strcmp(a->scID, b->scID)==0)
                {
                    ullong_save_add(qualTmp.sameUmiSameAbSc, 1);
                }
                else
                {
                    ullong_save_add(qualTmp.sameUmiDiffAbSc, 1);
                }
            }
        }
        ++tmpCurrentUmisChecked;
        if( (containerSize >= 100) && (tmpCurrentUmisChecked % (containerSize / 100) == 0) )
        {
            lock.lock();
            currentUmisChecked += tmpCurrentUmisChecked;
            tmpCurrentUmisChecked = 0;
            double perc = currentUmisChecked/ (double) containerSize;
            printProgress(perc);
            lock.unlock();
        }
    }
}

void BarcodeProcessingHandler::removeFalseSingleCellsFromUmis(const std::vector< std::vector<dataLinePtr> >& uniqueUmis, int& currentUmisChecked,
                                                                       std::vector<dataLinePtr>& dataLinesToDelete)
{
    // go over UMi data: safe all ABSc combos that have to be deleted
    //delete by Ab_id, sc_id and umi_id
    int containerSize = rawData.getUniqueUmis().size();
    int tmpCurrentUmisChecked =0;
    for(auto uniqueUmi : uniqueUmis)
    {
        //for all AbSc combinations of this unique UMI
        std::unordered_map<std::string, long> umiCountMap;
        for(int i = 0; i < uniqueUmi.size(); ++i)
        {
            std::unordered_map<std::string, long>::iterator umiCountMapIt = umiCountMap.find(uniqueUmi.at(i)->scID);
            if( umiCountMapIt != umiCountMap.end())
            {
                ++(umiCountMapIt->second);
            }
            else
            {
                umiCountMap.insert(std::make_pair(uniqueUmi.at(i)->scID, 1));
            }
        }

        std::string realSingleCellID;
        bool realSingleCellExists = false;
        for(auto singleCellCountPair : umiCountMap)
        {
            double singleCellPerc = singleCellCountPair.second/umiCountMap.size();
            if( singleCellPerc >= 0.9 )
            {
                realSingleCellID = singleCellCountPair.first;
                realSingleCellExists = true;
            }
        }

        if(!realSingleCellExists)
        {
            //for(int i = 0; i < uniqueUmi.size(); ++i)
            //{
             //   {
              //      dataLinesToDelete.push_back(uniqueUmi.at(i));
               // }
            //}
            continue;
        }
        else
        {
            for(int i = 0; i < uniqueUmi.size(); ++i)
            {
                if(uniqueUmi.at(i)->scID != realSingleCellID)
                {
                    dataLinesToDelete.push_back(uniqueUmi.at(i));
                }
            }

        }

        ++tmpCurrentUmisChecked;
        if( (containerSize >= 100) && (tmpCurrentUmisChecked % (containerSize / 100) == 0) )
        {
            lock.lock();
            currentUmisChecked += tmpCurrentUmisChecked;
            tmpCurrentUmisChecked = 0;
            double perc = currentUmisChecked/ (double) containerSize;
            printProgress(perc);
            lock.unlock();
        }
    }
}

void BarcodeProcessingHandler::writeStats(std::string output)
{
    //WRITE INTO FILE
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");
    if(found == std::string::npos)
    {
        output = "STATS" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "STATS" + output.substr(found+1);
    }
    outputFile.open (output);

    //qualiry output
    outputFile << "same umi same AbSc: " << qual.sameUmiSameAbSc << "\n";
    outputFile << "same umi diff AbSc: " << qual.sameUmiDiffAbSc << "\n";    

    //mismatches dist output
    for(auto el : stats.umiMismatchDict)
    {
        outputFile << el.first << " " << el.second << "\n";
    }

    outputFile.close();

}

void BarcodeProcessingHandler::writeUmiCorrectedData(const std::string& output)
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
    outputFile << "UMI" << "\t" << "AB_BARCODE" << "\t" << "SingleCell_BARCODE" << "\n"; 
    for(auto line : umiData)
    {
        outputFile << line->umiSeq << "\t" << line->abName << "\t" << line->scID << "\n"; 
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
    outputFile << "AB_BARCODE" << "\t" << "SingleCell_BARCODE" << "\t" << "AB_COUNT" << "\t" << "TREATMENT" << "\n"; 
    for(auto line : abData)
    {
        outputFile << (line.abName) << "\t" << line.scID << "\t" << line.abCount << "\t" << (line.treatment) << "\n"; 
    }
    outputFile.close();
}