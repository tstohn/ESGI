#include "UmiDataParser.hpp"
#include "helper.hpp"

void UmiDataParser::parseFile(const std::string fileName, const int& thread)
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

void UmiDataParser::parseBarcodeLines(std::istream* instream, const int& totalReads, int& currentReads)
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
            getCiBarcodeInWholeSequence(line, elements);
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

void UmiDataParser::addFastqReadToUmiData(const std::string& line, const int& elements)
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
    std::string treatment = "";

    if(treatmentIdx != INT_MAX)
    {
        treatment = result.at(treatmentIdx);
    }
    rawData.add(result.at(umiIdx), result.at(abIdx), singleCellIdx, treatment);

}

std::string UmiDataParser::generateSingleCellIndexFromBarcodes(std::vector<std::string> ciBarcodes)
{
    std::string scIdx;

    for(int i = 0; i < ciBarcodes.size(); ++i)
    {
        std::string barcodeAlternative = ciBarcodes.at(i);
        int tmpIdx = barcodeDict.barcodeIdDict.at(i)[barcodeAlternative];
        scIdx += std::to_string(tmpIdx);
        if(i < ciBarcodes.size() - 1)
        {
            scIdx += ".";
        }
    }

    return scIdx;
}

void UmiDataParser::getCiBarcodeInWholeSequence(const std::string& line, int& barcodeElements)
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
            if(variableBarcodeCount == barcodeDict.tmpTreatmentIdx)
            {
                treatmentIdx = count;
            }
            //add index for combinatorial indexing barcode
            if (std::count(barcodeDict.ciBarcodeIndices.begin(), barcodeDict.ciBarcodeIndices.end(), variableBarcodeCount)) 
            {
                fastqReadBarcodeIdx.push_back(count);
            }
            else
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
    assert(fastqReadBarcodeIdx.size() == barcodeDict.ciBarcodeIndices.size());
    assert(umiIdx != INT_MAX);
    assert(abIdx != INT_MAX);
}

void UmiDataParser::processBarcodeMapping(const int& umiMismatches, const int& thread)
{
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
    std::vector<std::vector<abLine> > abDataThreaded(thread);

    //for every batch calcualte the unique umis with same/different AbSc barcodes
    std::vector<std::thread> workers;
    int currentUmisChecked = 0;
    std::cout << "Checking Quality of UMIs\n";
    for (int i = 0; i < thread; ++i) 
    {
        workers.push_back(std::thread(&UmiDataParser::umiQualityCheck, this, std::ref(independantUmiBatches.at(i)), std::ref(umiQualThreaded.at(i)), std::ref(currentUmisChecked) ));
    }
    printProgress(1);
    for (std::thread &t: workers) 
    {
        if (t.joinable()) {
            t.join();
        }
    }
    std::cout << "\n";
    //for every batch calculate an UmiData and ABData vector and stats
    int currentUmisCorrected = 0;
    std::cout << "Correcting UMIs and counting ABs\n";
    for (int i = 0; i < thread; ++i) 
    {
        workers.push_back(std::thread(&UmiDataParser::correctUmis, this, std::ref(umiMismatches), std::ref(umiStatsThreaded.at(i)), std::ref(umiDataThreaded.at(i)),
                          std::ref(abDataThreaded.at(i)), std::ref(independantAbScBatches.at(i)), std::ref(currentUmisCorrected) ));
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

//correct for mismatches in the UMI
void UmiDataParser::correctUmis(const int& umiMismatches, StatsUmi& statsTmp, std::vector<dataLinePtr>& umiDataTmp, std::vector<abLine>& abDataTmp, 
                                const std::vector<std::vector<dataLinePtr> >& AbScBucket, int& currentUmisCorrected)
{
    //correct for UMI mismatches and fill the AbCountvector
    //iterate through same AbScIdx, calculate levenshtein dist for all UMIs and match those with a certain number of mismatches
    int tmpCurrentUmisCorrected = 0;
    int containerSize = rawData.getUniqueAbSc().size();
    for (auto uniqueAbSc : AbScBucket)
    {
        abLine abLineTmp;
        abLineTmp.cell_seq = uniqueAbSc.at(0)->cell_seq;
        abLineTmp.ab_seq = rawData.getProteinName(uniqueAbSc.at(0)->ab_seq);
        abLineTmp.treatment = rawData.getTreatmentName(uniqueAbSc.at(0)->treatment_seq);

        int abCount = 0; // abCount is calculated for every umi one after the other, if an umi is unique, the count is incremented
        //for umis with several occurences the last occurence will increment the count, the very last UMI is never checked and is always
        //incrementing the count, therefore initialized with 1
        
        while(!uniqueAbSc.empty())
        {
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
                const char* umia = uniqueAbSc.at(i)->umi_seq;
                const char* umib = uniqueAbSc.at(j)->umi_seq;
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
                            uniqueAbSc.at(j)->umi_seq = uniqueAbSc.at(i)->umi_seq;
                            rawData.changeUmi(umib, umia, uniqueAbSc.at(i));
                        }
                        else if(strlen(umib)==umiLength)
                        {
                            realUmi = umib;
                            uniqueAbSc.at(i)->umi_seq = uniqueAbSc.at(j)->umi_seq;
                            rawData.changeUmi(umia, umib, uniqueAbSc.at(j));
                        }
                        else
                        {
                            realUmi = umia;
                            uniqueAbSc.at(j)->umi_seq = uniqueAbSc.at(i)->umi_seq;
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

        abLineTmp.ab_cout = abCount;
        abDataTmp.push_back(abLineTmp);

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
void UmiDataParser::correctUmisWithStats(const int& umiMismatches, StatsUmi& statsTmp, std::vector<dataLinePtr>& umiDataTmp, std::vector<abLine>& abDataTmp, 
                                const std::vector<std::vector<dataLinePtr> >& AbScBucket, int& currentUmisCorrected)
{
    //correct for UMI mismatches and fill the AbCountvector
    //iterate through same AbScIdx, calculate levenshtein dist for all UMIs and match those with a certain number of mismatches
    int tmpCurrentUmisCorrected = 0;
    int containerSize = rawData.getUniqueAbSc().size();
    for (auto uniqueAbSc : AbScBucket)
    {
        abLine abLineTmp;
        abLineTmp.cell_seq = uniqueAbSc.at(0)->cell_seq;
        abLineTmp.ab_seq = rawData.getProteinName(uniqueAbSc.at(0)->ab_seq);
        abLineTmp.treatment = rawData.getTreatmentName(uniqueAbSc.at(0)->treatment_seq);

        int abCount = 1; // abCount is calculated for every umi one after the other, if an umi is unique, the count is incremented
        //for umis with several occurences the last occurence will increment the count, the very last UMI is never checked and is always
        //incrementing the count, therefore initialized with 1
        for(int i = 0; i < (uniqueAbSc.size() - 1); ++i)
        {
            //assert(abLineTmp.ab_seq == uniqueAbSc.second.at(i)->ab_seq);
            //assert(abLineTmp.cell_seq == uniqueAbSc.second.at(i)->cell_seq);
            bool unique = true;
            for(int j = i+1; j < uniqueAbSc.size(); ++j)
            {
                const char* umia = uniqueAbSc.at(i)->umi_seq;
                const char* umib = uniqueAbSc.at(j)->umi_seq;

                int dist = INT_MAX;
                int start = 0;
                int end = 0;

                bool similar = levenshtein(umia, umib, umiLength, start, end, dist, true);

                //if(std::strcmp(abLineTmp.cell_seq, "10295") == 0 & std::strcmp((abLineTmp.ab_seq)->c_str(), "CTD1")==0){std::cout << umia << " " << umib << " => "<< dist << " " << start << " " << end <<  "\n";}

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
                            uniqueAbSc.at(j)->umi_seq = uniqueAbSc.at(i)->umi_seq;
                            rawData.changeUmi(umib, umia, uniqueAbSc.at(i));
                        }
                        else if(strlen(umib)==umiLength)
                        {
                            realUmi = umib;
                            uniqueAbSc.at(i)->umi_seq = uniqueAbSc.at(j)->umi_seq;
                            rawData.changeUmi(umia, umib, uniqueAbSc.at(j));
                        }
                        else
                        {
                            realUmi = umia;
                            uniqueAbSc.at(j)->umi_seq = uniqueAbSc.at(i)->umi_seq;
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
            //if(std::strcmp(abLineTmp.cell_seq, "10295") == 0 & std::strcmp((abLineTmp.ab_seq)->c_str(), "CTD1")==0){std::cout << "\n";}
            if(unique)
            {
                ++abCount;
                            //if(std::strcmp(abLineTmp.cell_seq, "10295") == 0 & std::strcmp((abLineTmp.ab_seq)->c_str(), "CTD1")==0){std::cout << "      "<< abCount << "      ADDED\n";;}
            }
            umiDataTmp.push_back(uniqueAbSc.at(i));
        }    
        abLineTmp.ab_cout = abCount;
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

void UmiDataParser::umiQualityCheck(const std::vector< std::vector<dataLinePtr> >& uniqueUmis, umiQuality& qualTmp, int& currentUmisChecked)
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
                if(strcmp(a->ab_seq, b->ab_seq)==0 && strcmp(a->cell_seq, b->cell_seq)==0)
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

void UmiDataParser::writeStats(std::string output)
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

void UmiDataParser::writeUmiCorrectedData(const std::string& output)
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
        outputFile << line->umi_seq << "\t" << line->ab_seq << "\t" << line->cell_seq << "\n"; 
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
        outputFile << *(line.ab_seq) << "\t" << line.cell_seq << "\t" << line.ab_cout << "\t" << *(line.treatment) << "\n"; 
    }
    outputFile.close();
}