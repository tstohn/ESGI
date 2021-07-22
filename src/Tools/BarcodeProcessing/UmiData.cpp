#include "UmiData.hpp"
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
    
    //Cleanup
    file.close();
}

void UmiDataParser::parseBarcodeLines(std::istream* instream, const int& totalReads, int& currentReads)
{
    std::string line;
    std::cout << "READING ALL LINES INTO MEMORY\n";
    while(std::getline(*instream, line))
    {
        //for the first read check the positions in the string that refer to CIBarcoding positions
        if(currentReads==0){
            ++currentReads; 
            getCiBarcodeInWholeSequence(line);
            continue;
        }
        addFastqReadToUmiData(line);   

        double perc = currentReads/ (double)totalReads;
        ++currentReads;
        printProgress(perc);        
    }
    std::cout << "\n";
}

void UmiDataParser::addFastqReadToUmiData(const std::string& line)
{
    dataLine dataTmp;

    //split the line into barcodes
    std::vector<std::string> result;
    std::stringstream ss;
    ss.str(line);
    while( ss.good() )
    {
        std::string substr;
        getline( ss, substr, '\t' );
        result.push_back( substr );
    }

    //hand over the UMI string, ab string, singleCellstring (concatenation of CIbarcodes)
    std::vector<std::string> ciBarcodes;
    for(int i : fastqReadBarcodeIdx)
    {
        ciBarcodes.push_back(result.at(i));
    }
    std::string singleCellIdx = generateSingleCellIndexFromBarcodes(ciBarcodes);
    data.add(result.at(umiIdx), result.at(abIdx), singleCellIdx);

}

std::string UmiDataParser::generateSingleCellIndexFromBarcodes(std::vector<std::string> ciBarcodes)
{
    std::string scIdx;

    for(int i = 0; i < ciBarcodes.size(); ++i)
    {
        std::string barcodeAlternative = ciBarcodes.at(i);
        int tmpIdx = barcodeDict.barcodeIdDict.at(i)[barcodeAlternative];
        scIdx += std::to_string(tmpIdx);
    }

    return scIdx;
}

void UmiDataParser::getCiBarcodeInWholeSequence(const std::string& line)
{
    std::vector<std::string> result;
    std::stringstream ss;
    ss.str(line);
    int count = 0;
    int variableBarcodeCount = 0;
    while( ss.good() )
    {
        std::string substr;
        getline(ss, substr, '\t' );
        if(substr.empty()){continue;}
        //if substr is only N's
        if(substr.find_first_not_of('N') == std::string::npos)
        {
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
    assert(sizeof(fastqReadBarcodeIdx) == sizeof(barcodeDict.ciBarcodeIndices));
    assert(umiIdx != INT_MAX);
    assert(abIdx != INT_MAX);
}

void UmiDataParser::correctUmis(const int& umiMismatches)
{
    //first quality check, does a unique umi have always the same AbScIdx
    for(auto uniqueUmi : data.getUniqueUmis())
    {
        //for all AbSc combinations of this unique UMI
        for(int i = 0; i < (uniqueUmi.second.size() - 1); ++i)
        {
            for(int j = i+1; j < uniqueUmi.second.size(); ++j)
            {
                dataLinePtr a = uniqueUmi.second.at(i);
                dataLinePtr b = uniqueUmi.second.at(j);
                if(strcmp(a->ab_seq, b->ab_seq)==0 & strcmp(a->cell_seq, b->cell_seq)==0)
                {
                    ++stats.sameUmiSameAbSc;
                }
                else
                {
                    ++stats.sameUmiDiffAbSc;
                }
            }
        }
    }

    //correct for UMI mismatches
    //iterate through same AbScIdx, calculate levenshtein dist for all UMIs and match those with a certain number of mismatches
    for (auto uniqueAbSc : data.getUniqueAbSc())
    {
        for(int i = 0; i < (uniqueAbSc.second.size() - 1); ++i)
        {
            for(int j = i+1; j < uniqueAbSc.second.size(); ++j)
            {
                const char* umia = uniqueAbSc.second.at(i)->umi_seq;
                const char* umib = uniqueAbSc.second.at(j)->umi_seq;

                int dist = UINT_MAX;
                int start = 0;
                int end = 0;
                bool similar = levenshtein(umia, umib, umiLength, start, end, dist);

                //if mismatches are within range, change UMI seq
                //the new 'correct' UMI sequence is the one of umiLength, if both r of
                //same length, its the first occuring UMI
                if(dist <= umiMismatches & dist!=0)
                {
                    std::cout << "TWO SIMILAR ONES\n";
                    //get real UMI
                    const char* realUmi;
                    if(strlen(umia)==umiLength)
                    {
                        realUmi = umia;
                        uniqueAbSc.second.at(j)->umi_seq = uniqueAbSc.second.at(i)->umi_seq;
                        data.changeUmi(umib, umia, uniqueAbSc.second.at(i));
                    }
                    else if(strlen(umib)==umiLength)
                    {
                        realUmi = umib;
                        uniqueAbSc.second.at(i)->umi_seq = uniqueAbSc.second.at(j)->umi_seq;
                        data.changeUmi(umia, umib, uniqueAbSc.second.at(j));
                    }
                    else
                    {
                        realUmi = umia;
                        uniqueAbSc.second.at(j)->umi_seq = uniqueAbSc.second.at(i)->umi_seq;
                        data.changeUmi(umib, umia, uniqueAbSc.second.at(i));
                    }                    
                }

                //dist is set inside levenshtein, if its <= mismatches=2
                if(stats.umiMismatchDict.find(dist) != stats.umiMismatchDict.end())
                {
                    ++stats.umiMismatchDict.at(dist);
                }
                else
                {
                    stats.umiMismatchDict.insert(std::make_pair(dist, 1));
                }
            }
        }    
    }

}

void UmiDataParser::writeStats(std::string output)
{
    for(auto uniqueUmiIdx : data.getUniqueUmis())
    {
        if(uniqueUmiIdx.second.size() != 1)
        {
            std::cout << uniqueUmiIdx.first << " : " << uniqueUmiIdx.second.size() << "\n";
            auto duplicatedUmis = data.getDataWithUmi(uniqueUmiIdx.first);
            std::cout << "=>  ";
            for(auto umi : duplicatedUmis)
            {
                std::cout << umi->ab_seq << " " << umi->cell_seq << "; ";
            }
            std::cout << "\n";
        }
    }
    std::cout << "###################\n";
    for(auto uniqueAbScIdx : data.getUniqueAbSc())
    {
        if(uniqueAbScIdx.second.size() != 1)
        {
            std::cout << uniqueAbScIdx.first << " : " << uniqueAbScIdx.second.size() << "\n";
            auto duplicatedAbScs = data.getDataWithAbSc(uniqueAbScIdx.first);
            std::cout << "=>  ";
            for(auto absc : duplicatedAbScs)
            {
                std::cout << absc->umi_seq << " " << absc->ab_seq << " " << absc->cell_seq << "; ";
            }
            std::cout << "\n";
        }
    }

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
    outputFile << "same umi same AbSc: " << stats.sameUmiSameAbSc << "\n";
    outputFile << "same umi diff AbSc: " << stats.sameUmiDiffAbSc << "\n";    

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
    outputFile.open (output);

    outputFile << "UMI" << "\t" << "AB_BARCODE" << "\t" << "SingleCell_BARCODE" << "\n"; 
    for(auto line : data.getData())
    {
        outputFile << line->umi_seq << "\t" << line->ab_seq << "\t" << line->cell_seq << "\n"; 
    }

    outputFile.close();
}