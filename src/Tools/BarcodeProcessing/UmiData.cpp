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
    while(std::getline(*instream, line))
    {
        if(currentReads==0){
            ++currentReads; 
            getCiBarcodeInWholeSequence(line);
            continue;
        }
        std::cout << "adding line " <<currentReads <<  line << "\n";
        addFastqReadToUmiData(line);   
        std::cout << "SUCCESS for " << line << "\n";

        double perc = currentReads/ (double)totalReads;
        ++currentReads;
        //printProgress(perc);        
    }
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
    std::cout << "make SCINDEX\n";
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
                std::cout << "barcode " << count << "\n";
            }
            else
            {
                abIdx = count;
            std::cout << "AB"<<abIdx << "\n";

            }
            ++variableBarcodeCount;
        }
        else if(substr.find_first_not_of('X') == std::string::npos)
        {
            umiIdx = count;
            std::cout << "UMI"<< umiIdx << "\n";
        }
        result.push_back( substr );
        ++count;
    }
    assert(sizeof(fastqReadBarcodeIdx) == sizeof(barcodeDict.ciBarcodeIndices));
    assert(umiIdx != INT_MAX);
    assert(abIdx != INT_MAX);
}