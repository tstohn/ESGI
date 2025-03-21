#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <string_view>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <unordered_map>
#include <mutex>

struct fastqLine
{
    std::string line;
    std::string quality = "";
    std::string name = "";
};

class DemultiplexedLine
{
    public:

        std::vector<std::string> barcodeList;

        bool containsDNA = false;
        std::string dna = "";
        std::string dnaQuality = "";
        std::string dnaName = "";
};

typedef std::vector<const char*> BarcodeMapping;
typedef std::vector<BarcodeMapping> BarcodeMappingVector;
/** @brief representation of all the mapped barcodes:
 * basically a vector of all reads, where each read itself is a vector of all mapped barcodes
 * This structures stores each barcode only once, handled by the UniqueCharSet, by that
 * most highly redundant datasets can be stored in only a fraction of its origional memory
**/
class DemultiplexedReads
{
    public:

        DemultiplexedReads()
        {
            uniqueChars = std::make_shared<UniqueCharSet>();
            lock = std::make_unique<std::mutex>();
        }

        void addVector(std::vector<std::string> barcodeVector)
        {
            std::lock_guard<std::mutex> guard(*lock);
            BarcodeMapping uniqueBarcodeVector;
            for(std::string barcode : barcodeVector)
            {
                uniqueBarcodeVector.emplace_back(uniqueChars->getUniqueChar(barcode.c_str()));
            }
            mappedBarcodes.push_back(uniqueBarcodeVector);
        }

        const size_t size()
        {
            return mappedBarcodes.size();
        }

        const BarcodeMapping at(const int& i)
        {
            return(mappedBarcodes.at(i));
        }

        const BarcodeMappingVector get_all_reads()
        {
            return(mappedBarcodes);
        }

    private:
        BarcodeMappingVector mappedBarcodes;
        //all the string inside this class are stored only once, 
        //set of all the unique barcodes we use, and we only pass pointers to those
        std::shared_ptr<UniqueCharSet> uniqueChars;
        std::unique_ptr<std::mutex> lock;  

};

typedef std::shared_ptr<DemultiplexedReads> DemultiplexedReadsPtr; 
