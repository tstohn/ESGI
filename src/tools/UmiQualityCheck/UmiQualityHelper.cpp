#include "UmiQualityHelper.hpp"


void UmiQuality::writeUmiQualityData(std::string output)
{
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");

    //STORE RAW UMI CORRECTED DATA
    std::string umiOutput = output;
    if(found == std::string::npos)
    {
        umiOutput = "UmiQualityCheck" + output;
    }
    else
    {
        umiOutput = output.substr(0,found) + "/" + "UmiQualityCheck" + output.substr(found+1);
    }
    outputFile.open (umiOutput);

    int count = 0;
    for(auto key : umiQualStat.get_keys_in_order())
    {
        outputFile << key;
        if(count==umiQualStat.get_keys_in_order().size()-1)
        {
            outputFile << "\t" << "COUNT\n";
        }
        else
        {
            outputFile << "\t";
        }
         ++count;
    }
    for(auto errorMapElem : umiQualStat.get_error_map())
    {
        count = 0;
        for(auto barcodeCounts : errorMapElem.first)
        {
            if(count==errorMapElem.first.size()-1)
            {
                outputFile << barcodeCounts << "\t" << errorMapElem.second << "\n";
            }
            else
            {
                outputFile << barcodeCounts << "\t";
            }
            ++count;
        }

    }
    outputFile.close();
}

void addToMap(const std::string& key, const std::string& value, std::unordered_map<std::string, std::vector<std::string> >& uniqueBarcodes)
{
    std::unordered_map<std::string, std::vector<std::string> >::iterator uniqueBarcodesIt = uniqueBarcodes.find(key);
    //if the key already exists add to it
    if( uniqueBarcodesIt != uniqueBarcodes.end())
    {
        //add value only if it is not already inside
        if (!std::count(uniqueBarcodesIt->second.begin(), uniqueBarcodesIt->second.end(), value)) 
        {
            uniqueBarcodesIt->second.push_back(value);
        }
    }
    else
    {
        std::vector<std::string> valueVec;
        valueVec.push_back(value);
        uniqueBarcodes.insert(std::make_pair(key, valueVec));
    }
}

void umiQualityStat::add_value(const std::unordered_map<std::string, std::vector<std::string> >& barcodesForUmi)
{
    //generate key (number of different barcodes for each barcode(AB, BC1, BC2, ...) )
    std::unordered_map<std::string, int> countMap; //map storing how many different barcodes we see in each BC-round(Ab,BC1, BC2,...)
    for(const std::pair<std::string, std::vector<std::string>>& elem : barcodesForUmi)
    {
        int count = static_cast<int>(elem.second.size());
        countMap.insert(std::make_pair(elem.first, count));
    }

    std::lock_guard<std::mutex> guard(errorMapUpdateLock);
    
    if(errorMap.empty())
    {
        //add order of keys to 'barcodeOrder'
        for(auto elem : countMap) 
        {
            barcodeOrder.push_back(elem.first);
        }
    }

    std::vector<int> barcodeCountVector;
    for(auto barcodestring : barcodeOrder)
    {
        barcodeCountVector.push_back(countMap.at(barcodestring));
    }
    std::unordered_map< std::vector<int>, unsigned long long, VectorHasher>::iterator errorMapIt = errorMap.find(barcodeCountVector);
    if( errorMapIt != errorMap.end())
    {
        ullong_save_add(errorMapIt->second, 1);
    }
    else
    {
        errorMap.insert(std::make_pair(barcodeCountVector, 1));
    }
}


void UmiQuality::checkUniquenessOfUmis(const std::vector<dataLinePtr>& uniqueUmiLines)
{
    std::unordered_map<std::string, std::vector<std::string> > uniqueBarcodes; //maps string for Bc type (AB,BC1,...) => vector of all the possible barcode that we encounter

    //analyze each dataLine and store unique occurences of each barcode
    for(const dataLinePtr& line : uniqueUmiLines)
    {
        //for all CI barcodes
        std::vector<std::string> ciVec = splitByDelimiter(line->scID, ".");
        int CiBcCount = 0;
        for(const std::string& ciBarcodeId : ciVec)
        {
            std::string key = "CiIdx" + std::to_string(CiBcCount);
            addToMap(key, ciBarcodeId, uniqueBarcodes);
            ++CiBcCount;
        }

        //for Ab barcode
        std::string abKey = "AB";
        addToMap(abKey, line->abName, uniqueBarcodes);

        //for Treatment barcode
        std::string treatKey = "TREATMENT";
        addToMap(treatKey, line->treatmentName, uniqueBarcodes);
    }

    //then safely add this information to 'umiQualityStat'
    umiQualStat.add_value(uniqueBarcodes);
}

void UmiQuality::runUmiQualityCheck(const int& thread, const std::string& output)
{
    boost::asio::thread_pool pool(thread); //create thread pool
    //Map of UMIs with duplicate Sc-Ab-treatments
    for(std::pair<const char *, std::vector<dataLinePtr>> dataLinesOfUmi : *rawData.getUniqueUmis())
    {
        //handing over only lineCount as reference, everything else will be copied (Mapping object as handed overr as this-pointer)
        boost::asio::post(pool, std::bind(&UmiQuality::checkUniquenessOfUmis, this, dataLinesOfUmi.second));
    }
    pool.join();

    writeUmiQualityData(output);
}