#include "BarcodeProcessingHandler.hpp"

#include <unordered_set>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <filesystem>

std::string stripExtension(const std::string& filename) 
{
    std::filesystem::path p(filename);
    if (p.has_extension()) {
        return p.stem().string();  // removes last extension
    }
    return filename;  // no extension, return unchanged
}

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
    if (!barcodeFileStream.is_open()) {
        std::cerr << "Error: Could not open following file for barcode parsing: " << file << std::endl;
        std::cerr << "Please double check if the path to the barcode-files is right." << std::endl;
        // check wrong spaces
        if (file.find(',') != std::string::npos || file.find(' ') != std::string::npos || file.find('\n') != std::string::npos) 
        {
            std::cout << "The file path contains spacing characters like: <,> or < > or <newline>. \n";
            std::cout << "Did you accidentally use one of them instead of tabs to seperate columns.\n";
        }
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::vector<std::string> seqVector;
    while (std::getline(barcodeFileStream, line)) 
    {
        //check Windows-specific trailing newlines
        if (!line.empty() && (line.back() == '\n' || line.back() == '\r')) 
        {
            line.pop_back();
        }

        std::string delimiter = ",";
        std::string seq;
        size_t pos = 0;
        while ((pos = line.find(delimiter)) != std::string::npos) 
        {
            seq = trim(line.substr(0, pos));
            line.erase(0, pos + 1);
            for (char const &c: seq) {
                if(!( (c=='A') | (c=='T') | (c=='G') | (c=='C') |
                        (c=='a') | (c=='t') | (c=='g') | (c=='c')))
                        {
                        std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C)\n";
                        if( (c==' ') | (c=='\t') | (c=='\n'))
                        {
                            std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
                        }
                        exit(1);
                        }
            }
            seqVector.push_back(seq);
        }
        seq = trim(line);
        for (char const &c: seq) {
            if(!((c=='A') || (c=='T') || (c=='G') || (c=='C') ||
                    (c=='a') || (c=='t') || (c=='g') || (c=='c')))
                    {
                    std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C)\n";
                    if((c==' ') || (c=='\t') || (c=='\n'))
                    {
                        std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
                    }
                    exit(1);
                    }
        }
        if(seq!="")
        {
            seqVector.push_back(seq);
        }
    }

    barcodeList.insert(std::make_pair(colIdx, seqVector));
    seqVector.clear();
}

std::vector<std::string> parseFeatureNames(const std::string& featureFile) 
{
    std::vector<std::string> result;
    std::stringstream ss(featureFile);
    std::string item;
    while (std::getline(ss, item, ',')) 
    {
        //check Windows-specific trailing newlines
        if (!item.empty() && (item.back() == '\n' || item.back() == '\r')) 
        {
            item.pop_back();
        }
        result.push_back(item);
    }
    return result;
}

void generateBarcodeDicts(const std::string& headerLine, const std::string& barcodeDir, std::string barcodeIndices, BarcodeInformation& barcodeIdData, 
                          std::vector<std::string>& proteinNamelist, bool parseAbBarcodes, const int& featureIdx, bool& umiRemoval,
                          std::vector<std::vector<std::string>>* annotationBarcodesVector, const std::vector<int>* annotationIdxs, std::string umiIdx,
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
    std::unordered_map<int, std::vector<std::string>> barcodeList; //maps column number -> all barcodes: e.g. 2 -> variables barcodes in BC1.txt, ...
    while (std::getline(ss, colName, '\t')) 
    {
        //check Windows-specific trailing newlines
        if (!colName.empty() && (colName.back() == '\n' || colName.back() == '\r')) 
        {
            colName.pop_back();
        }

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

        //test validity of single-cell indices (do we have enough headers)
        int max_scIdx = barcodeIdData.scBarcodeIndices.empty() 
                            ? 0 
                            : *std::max_element(barcodeIdData.scBarcodeIndices.begin(), 
                                                barcodeIdData.scBarcodeIndices.end());
        if (max_scIdx >= colIdx)
        {
            std::cout << "Double check the input file and indices for single-cells.\n";
            std::cout << "There is at least one single-cell column that does not exist (higher single-cell index than number of headers)!\n";
            std::cout << "Be aware that the first column might contain the readname, and columns are 0-indexed!\n";
            exit(EXIT_FAILURE);
        }

        std::cout << "Assigning single-cells according to columns: ";
        for(size_t i = 0; i < barcodeIdData.scBarcodeIndices.size()-1; ++i)
        {
            std::cout << barcodeHeader.at(barcodeIdData.scBarcodeIndices.at(i)) << ",";

            if(!endsWithTxt(barcodeHeader.at(barcodeIdData.scBarcodeIndices.at(i))))
            {
                std::cout << "\nError when parsing the potential barcodes for single-cell column <" << barcodeHeader.at(barcodeIdData.scBarcodeIndices.at(i)) << ">\n";
                std::cout << "Unfortunatelly, the file with single-cell barcodes MUST end with .txt and should contain comma seperated barcodes.\n";
                std::cout << "The supplied file seems to not end with .txt\n";

                exit(EXIT_FAILURE);
            }
        }
        std::cout << barcodeHeader.at(barcodeIdData.scBarcodeIndices.back());
        std::cout << "\n";
        if(!endsWithTxt(barcodeHeader.at(barcodeIdData.scBarcodeIndices.back())))
        {
            std::cout << "\nError when parsing the potential barcodes for single-cell column <" << barcodeHeader.at(barcodeIdData.scBarcodeIndices.back()) << ">\n";
            std::cout << "Unfortunatelly, the file with single-cell barcodes MUST end with .txt and should contain comma seperated barcodes.\n";
            std::cout << "The supplied file seems to not end with .txt\n";

            exit(EXIT_FAILURE);
        }

        //make map colIdx -> (map: barcode string -> number)
        //iterate over colIdx for scID
        for(unsigned int scIdx : barcodeIdData.scBarcodeIndices)
        {
            int barcodeCount = 0;
            std::unordered_map<std::string, int> barcodeMap;
            //map all barcode strings to numbers
            try
            {
                for(const std::string& barcodeEntry : barcodeList.at(scIdx))
                {
                    barcodeMap.insert(std::pair<std::string, int>(barcodeEntry, barcodeCount));
                    ++barcodeCount;
                }
            }catch (const std::out_of_range& e) 
            {
                std::cout << "Could not retrieve dictionaries of possible single-cell barcodes\n";
                std::cout << "Please ensure that the lists of single-cell barcodes can be correctly parsed\n";
                std::cout << "Therefore, make sure that the directory for barcode files exists (-d), \n";
                std::cout << "and that the files containing the barcodes are in this directory.\n";
                std::cout << "The file names must be the headers of the input demultiplexed tsv-file\n";
                exit(EXIT_FAILURE);
            }

            barcodeIdData.barcodeIdMaps.push_back(barcodeMap);
        }
    }

    //print assigned grouping index
    if(annotationIdxs != nullptr)
    {
        std::vector<int>::const_iterator it = std::max_element(annotationIdxs->begin(), annotationIdxs->end());
        if (it != annotationIdxs->end() && (*it) >= colIdx)
        {
            std::cout << "Double check the input file and annotation indices.\n";
            std::cout << "At least one annotation index is bigger than the number of headers in the input file!\n";
            exit(EXIT_FAILURE);
        }
        for(int aIdx : (*annotationIdxs))
        {
            std::cout << "Assigning annotations to column: " << barcodeHeader.at(aIdx) << "\n";
        }
        barcodeIdData.annotationIdxs = *annotationIdxs;  
        for(int aIdx : (*annotationIdxs))
        {
            std::string annotationString = "ANNOTATION_BC:" + stripExtension(barcodeHeader.at(aIdx)) + "_COL:" + std::to_string(aIdx);
            barcodeIdData.annotationFileCol.push_back(annotationString);
        }
    }

    //print the used feature index
    if (featureIdx >= colIdx)
    {
        std::cout << "Double check the input file and feature index.\n";
        std::cout << "The feature index is bigger than the number of headers in the input file!\n";
        exit(EXIT_FAILURE);
    }
    std::cout << "Assigning feature to column: " << barcodeHeader.at(featureIdx) << "\n";
    barcodeIdData.featureIdx = featureIdx;

    //print the used UMI indices
    if(barcodeIdData.umiIdx.empty())
    {
        std::cout << "Data contains NO UMI - we skip UMI-collapsing and count all reads!\n";
        umiRemoval = false;
        //umiRemoval is passed by reference, outside this function in main we will set the umiRemoval-parameter for the final counting-class - the BarcodeProcessingHandler
    }
    else if(umiRemoval == false)
    {
        std::cout << "You explicitely set umi-collpasing to false!\n";
        std::cout << "Skipping umi-collapsing deliberately and counting all reads\n";
    }
    else
    {
        std::cout << "Using following columns as UMIs: ";
        for(int colIdxForUMI : barcodeIdData.umiIdx)
        {
            std::cout << barcodeHeader.at(colIdxForUMI) << ",";
        }
        std::cout << "\n";
        std::cout << "Aligning all UMI-sequences with " << std::to_string(barcodeIdData.umiMismatches) << " mismatches.\n";
    }

    //assign the final dictionaries of variable barcodes
    if(parseAbBarcodes)
    {
        try 
        {
            proteinNamelist = barcodeList.at(featureIdx);
        } catch (const std::out_of_range& e) 
        {
            std::cout << "Please double check your input: Especially where is the feature column\n";
            throw std::runtime_error("There is no match for feature names in the assigned feature-barcode column: " + barcodeHeader.at(featureIdx));
        }
    }
    if(annotationIdxs != nullptr)
    {
        for(int aIdx : *annotationIdxs)
        {
            try 
            {
                annotationBarcodesVector->push_back(barcodeList.at(aIdx));
            } catch (const std::out_of_range& e) 
            {
                std::cout << "Please double check your input: Especially if the annotation columns are valid\n";
                throw std::runtime_error("There is no match for barcodes for the annotation-barcode column: " + barcodeHeader.at(aIdx));
            }
        }
    }

}

void BarcodeProcessingHandler::parse_barcode_sharing_file(std::string& barcodeFuseFile)
{
    std::ifstream in(barcodeFuseFile);
    if (!in) 
    {
        std::cerr << "Cannot open barcode fusing file: " << barcodeFuseFile << "\n";
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(in, line)) 
    {
        if (line.empty()) continue;

        std::istringstream stream(line);
        unsigned int position;
        std::string key, value;

        if (!(stream >> position >> value >> key)) 
        {
            std::cerr << "Skipping malformed line in barcode fuse file: " << line << "\n";
            continue;
        }

        // Insert or update the nested map
        barcodeSharingMap[position][key] = value;
    }
}

void BarcodeProcessingHandler::parse_barcode_file(const std::string& inFile)
{
    //reopen file for line counting
    std::ifstream file;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    bool gz = isGzipped(inFile);
    std::istream* instream = openFile(inFile, file, inbuf, gz);
    if (!instream)
    {
        std::cerr << "Error reading input file or file is empty! Please double check if the file exists:" << inFile << std::endl;
        exit(EXIT_FAILURE);
    }

    unsigned long long totalReads = 0;
    std::string line;
    while (std::getline(*instream, line)) 
    {
        ++totalReads;
    }
    if (instream != &file) delete instream;
    instream = nullptr;

    unsigned long long currentReads = 0;
    parseBarcodeLines(inFile, totalReads, currentReads);
}

void BarcodeProcessingHandler::parseBarcodeLines(const std::string& inFile, const unsigned long long& totalReads, unsigned long long& currentReads)
{
    //reopen file
    std::ifstream file;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    bool gz = isGzipped(inFile);
    std::istream* instream = openFile(inFile, file, inbuf, gz);
    if (!instream)
    {
        std::cerr << "Error reading input file or file is empty! Please double check if the file exists:" << inFile << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::cout << "STEP[1/3]\t(READING ALL LINES INTO MEMORY)\n";
    int elements = 0; //check that each row has the correct number of barcodes
    unsigned long long readCount = 0;
    while(std::getline(*instream, line))
    {
        //check Windows-specific trailing newlines
        if (!line.empty() && (line.back() == '\n' || line.back() == '\r')) 
        {
            line.pop_back();
        }

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

        ++currentReads;

        double perc = currentReads/ (double)totalReads;
        if(std::fmod(perc, 2.0) == 0.0)
        {
            printProgress(perc);        
        }
    }

    result.set_total_reads(currentReads-1); //minus header line
    result.set_total_ab_reads(readCount);

    printProgress(1);
    std::cout << "\n";
}

void BarcodeProcessingHandler::add_line_to_temporary_data(const std::string& line, const size_t& elements, unsigned long long& readCount)
{
    //split the line into barcodes
    std::vector<std::string> result;
    std::stringstream ss;
    ss.str(line);
    std::string substr;

    unsigned int position = 0;
    while(getline(ss, substr, '\t'))
    {
        if(substr != "")
        {
            //check if we have to replace barcodes
            if(!barcodeSharingMap.empty())
            {
                auto positionIt = barcodeSharingMap.find(position);
                //position has replacements
                if (positionIt != barcodeSharingMap.end()) 
                {
                    auto& barcodeMap = positionIt->second;
                    auto barcodeSubstrIt = barcodeMap.find(substr);
                    if (barcodeSubstrIt != barcodeMap.end()) 
                    {
                        substr = barcodeSubstrIt->second;
                    }
                }
            }
            result.push_back(substr);
            ++position;
        }
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
    
    std::vector<std::string> annotations;
    if(!barcodeInformation.annotationIdxs.empty())
    {
        for(int aIdx : barcodeInformation.annotationIdxs)
        {
            annotations.push_back(rawData.getAnnotationName(aIdx ,result.at(aIdx)));
        }
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

        rawData.add_to_umiDict(umiSeq, featureName, singleCellIdx, annotations);
    }
    //otherwise add reads directly to dict of ScAb to reads
    else
    {
        rawData.add_to_scAbDict("", featureName, singleCellIdx, annotations);
    }
}

std::string BarcodeProcessingHandler::generateSingleCellIndexFromBarcodes(const std::vector<std::string>& ciBarcodes)
{
    std::string scIdx;
    if(scIdString)
    {
        for(size_t i=0; i < ciBarcodes.size(); ++i)
        {
            scIdx += ciBarcodes.at(i);
        }
        return scIdx;
    }
    for(size_t i = 0; i < barcodeInformation.scBarcodeIndices.size(); ++i)
    {
        //first is the column idx of the barcode, second is the actual barcode that we want to map to a number
        int tmpIdx;
        //if we can not access this map, the barcode is invalid and has no index assigned
        try 
        {
            tmpIdx = (barcodeInformation.barcodeIdMaps.at(i)).at(ciBarcodes.at(i));
        } catch (const std::out_of_range& e) 
        {
            std::cout << "\nIt seems like there is an invalid barcode in the input table: " << ciBarcodes.at(i) << "\n";
            std::cout << "Please double check your input barcodes\n";
            throw std::runtime_error("Barcode not found in map assigning single-cell indices to barcodes: " + ciBarcodes.at(i));
            exit(EXIT_FAILURE);
        }

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
        if( singleCellPerc >= umiFilterThreshold) //default = 0.0
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
                                                   const std::vector<dataLinePtr>& allScAbCounts)
{
    //this stores only the number of UMIs that were collapsed into each other due to MM
    unsigned long long numberAlignedUmis = 0;
    //skip the first UMI, this is the one we compare all others to
    const char* umia = allScAbCounts.front()->umiSeq; //comapre first element to others
    for(size_t j = 1; j < allScAbCounts.size(); ++j)
    {
        //calling outputSense algorithm, much faster than levenshtein O(e*max(m,n))
        //however is recently implemented without backtracking
        //before umiMismatches was increased by the length difference between the two UMIs 
        //(no longer done, those deletion should probably be considered as part of the allowed umiMismatches)
        const char* umib = allScAbCounts.at(j)->umiSeq;
        unsigned int dist = UINT_MAX;
        bool similar;
        //the lower abundance UMI MUST BE represented less than 20% of the high abundance (first) UMI
        if( barcodeInformation.umiAbundanceThreshold > 0.0 && !(allScAbCounts.at(j)->umiCount < barcodeInformation.umiAbundanceThreshold * allScAbCounts.front()->umiCount))
        {
            similar = false;
        }
        else if(barcodeInformation.hamming)
        {
            similar = hamming_dist(umia, umib, barcodeInformation.umiMismatches, dist);
        }
        else
        {
            similar = outputSense(umia, umib, barcodeInformation.umiMismatches, dist);
        }

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

// collapse UMIs in vector of datalinePtrs, comapres only const char* by default
// make sure that umis Are stored as const char* where SAME Umis are stored in the same address
void BarcodeProcessingHandler::collapse_identical_UMIs(std::vector<dataLinePtr>& scAbCounts) 
{
    //iterate through UMIs and store positions that should be removed afterwards
    //DO NOT REMOVE current line and inside umiCount of the line store how many lines were collapsed
    for(size_t i = 0; i < (scAbCounts.size()-1); ++i)
    {
        std::vector<unsigned int> deletePositions;
        for(size_t j = i+1; j < scAbCounts.size(); ++j)
        {
            //check if UMIs are identical (were stored under the same char*)
            if(scAbCounts.at(i)->umiSeq == scAbCounts.at(j)->umiSeq)
            {
                //store positions to remove
                deletePositions.push_back(j);
                //increase UMI count of the 'master-line'
                ++(scAbCounts.at(i)->umiCount); // increase count for this UMI
            }
        }
        //remove all the positions that were collapsed and do not ahve to be checked anymore
        for(int posIdx = static_cast<int>(deletePositions.size()) - 1; posIdx >= 0; --posIdx)
        {
            unsigned int pos = deletePositions.at(posIdx);
            scAbCounts.erase(scAbCounts.begin() + pos);
        }
    }
}

void BarcodeProcessingHandler::count_abs_per_single_cell(const std::vector<dataLinePtr>& uniqueAbSc,
                                                        std::atomic<unsigned long long>& count,
                                                        const unsigned long long& totalCount)
{
        //correct for UMI mismatches and fill the AbCountvector
        //iterate through same AbScIdx, calculate levenshtein dist for all UMIs and match those with a certain number of mismatches

        //all dataLines for this AB SC combination
        std::vector<dataLinePtr> scAbCounts = uniqueAbSc;

        //data structures to be filled for the UMI and AB count
        scAbCount abLineTmp; // we fill only this one AB SC count
        umiCount umiLineTmp; //when iterating through scAbCounts the UMI is set new every time we encouter a new UMI

        abLineTmp.scID = umiLineTmp.scID = uniqueAbSc.at(0)->scID;
        abLineTmp.abName = umiLineTmp.abName = uniqueAbSc.at(0)->abName;
        abLineTmp.annotations = umiLineTmp.annotations = uniqueAbSc.at(0)->annotations;
        abLineTmp.className = uniqueAbSc.at(0)->cellClassname;

        //if we have no umis erase whole vector and count every element
        if(umiRemoval==false)
        {
            abLineTmp.abCount = scAbCounts.size();
            scAbCounts.clear();
        }
        else
        {
            //erase dataLinePtrs for lines with EXACTLY the same UMI (no MM) and increase umi count in this linePtr
            //can be done fast due to simple const char* comparison (we stored unique UMIs previously)
            collapse_identical_UMIs(scAbCounts);

            //new sort the remaining UMIs (after collpasing EXACTLY UNIQUE ones) in order of occurences, 
            //then start with the first (MOST ABUNDANT) for further demultiplexing
            std::sort(scAbCounts.begin(), scAbCounts.end(), sort_descending_by_umi_count);
        }

        //we take always first element in vector of read of same AB and SC ID (the element of most UMI counts)
        //then store all reads where UMIs are within distance, and delete those lines, and sum up their UMI counts (they might have been collapsed before on EXACT IDENTITY)
        while(!scAbCounts.empty())
        {
            dataLinePtr firstAbSc = scAbCounts.front();
            
            umiLineTmp.umi =  firstAbSc->umiSeq;
            umiLineTmp.abCount = firstAbSc->umiCount; //count the first occurence

            //check if we have to delete element anyways, bcs umi is too long
            //if( std::abs(int( strlen(lastAbSc->umiSeq) - barcodeInformation.umiLength )) >  barcodeInformation.umiMismatches)
            //{
                //scAbCounts.pop_back();
                //continue;
            //}

            //if in last element
        //    if(scAbCounts.size() == 1)
        //    {
         //       ++abLineTmp.abCount;
         //       umiLineTmp.abCount += firstAbSc->umiCount; 
         //       umiLineTmp.umi = firstAbSc->umiSeq;

        //        result.add_umi_count(umiLineTmp);
        //        result.add_umi_stats(umiLineTmp);
        //        scAbCounts.pop_back();
        //        break;
        //    }

            //otherwise conmpare all and mark the ones to delete
            std::vector<int> deletePositions;
            deletePositions.push_back(0); //add the first line to lines to delete

            //count all occurences of the first UMI for this AB-SC:
            // store the positions of other UMIs within umiMismatches to delete and not count those to final sc AB count (collapse them)
            if(barcodeInformation.umiMismatches > 0)
            {
                //add positions that should be deleted bcs. they contains same UMI, increase the count for this UMI in umiLineTmp
                count_umi_occurence(deletePositions, umiLineTmp, scAbCounts);
            }

            //ADD UMI if exists
            if(umiLineTmp.abCount > 0)
            {
                result.add_umi_count(umiLineTmp);
                result.add_umi_stats(umiLineTmp);
            }

            //delete all same UMIs
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

        if((totalCount >= 100) && ( ((100*count) / totalCount) % 5 == 0))
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
            
    const std::shared_ptr< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator>> AbScMap = rawData.getUniqueAbSc();
    for(std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator>::const_iterator it = AbScMap->begin(); 
        it != AbScMap->end(); 
        it++)
    {
        //as above: abSc is copied only
        boost::asio::post(pool_3, std::bind(&BarcodeProcessingHandler::count_abs_per_single_cell, this, 
                                            std::cref(it->second), std::ref(umiCount), std::cref(totalCount)));
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
        umiOutput = "UMIDATA_" + output;
    }
    else
    {
        umiOutput = output.substr(0,found) + "/" + "UMIDATA_" + output.substr(found+1);
    }
    outputFile.open (umiOutput);

    //write headers
    outputFile << "UMI" << "\t" << "FEATURE_ID" << "\t" << "SC_ID" << "\t";
    for(const std::string& annotationHeader : barcodeInformation.annotationFileCol)
    {
        outputFile << annotationHeader << "\t";
    }
    outputFile << "UMI_COUNT" << "\n"; 

    //write lines
    for(umiCount line : result.get_umi_data())
    {
        outputFile << line.umi << "\t" << line.abName << "\t" << line.scID << "\t";
        for(const char* annotation : line.annotations)
        {
            outputFile << annotation << "\t";
        }
        outputFile << line.abCount << "\n"; 
    }
    outputFile.close();

    //STORE AB COUNT DATA
    std::string abOutput = output;
    if(found == std::string::npos)
    {
        abOutput = "COUNTDATA_" + output;
    }
    else
    {
        abOutput = output.substr(0,found) + "/" + "COUNTDATA_" + output.substr(found+1);
    }
    outputFile.open (abOutput);
    bool writeClassLabels = rawData.check_class();
    if(writeClassLabels)
    {
        outputFile << "AB_BARCODE" << "\t" << "SC_BARCODE" << "\t" << "AB_COUNT" << "\t";
        for(const std::string& annotationHeader : barcodeInformation.annotationFileCol)
        {
            outputFile << annotationHeader << "\t";
        }
        outputFile  << "CLASS" << "\t" << "CLASS_COUNT" <<"\n"; 
    }
    else
    {
        outputFile << "FEATURE_ID" << "\t" << "SC_ID" << "\t" << "FEATURE_COUNT";
        //if there are no annotations the line ends here, otherwise continue with a tab...
        if(barcodeInformation.annotationFileCol.empty())
        {
            outputFile << "\n";
        }
        else
        {
            outputFile << "\t";
        }
        size_t count = 0;
        for(const std::string& annotationHeader : barcodeInformation.annotationFileCol)
        {
            ++count;
            outputFile << annotationHeader;
            if(count == barcodeInformation.annotationFileCol.size())
            {
                outputFile << "\n";
            }
            else
            {
                outputFile << "\t";
            }
        }
    }
    for(scAbCount line : result.get_ab_data())
    {
       /* if(writeClassLabels)
        {
                outputFile << line.abName << "\t" << line.scID << "\t" << line.abCount << "\t";
                 << line.treatment << "\t"
                 
                 << line.className << "\t" << guideCountPerSC.at(line.scID) << "\n"; 
        }*/

                outputFile << line.abName << "\t" << line.scID << "\t" << line.abCount;
                //write newline if we have no annotations
                if(line.annotations.empty())
                {
                    outputFile << "\n";
                }
                else
                {
                    outputFile << "\t";
                }
                size_t count = 0;
                for(const char* annotation : line.annotations)
                {
                    ++count;
                    outputFile << annotation; 
                    if(count == line.annotations.size())
                    {
                        outputFile << "\n";
                    }
                    else
                    {
                        outputFile << "\t";
                    }

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
