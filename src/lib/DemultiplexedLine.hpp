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

        std::string dna = "";
        std::string dnaQuality = "";
        std::string readName = "";
};

/** @brief representation of all the mapped barcodes:
 * basically a vector of all reads, where each read itself is a vector of all mapped barcodes
 * This structures stores each barcode only once, handled by the UniqueCharSet, by that
 * most highly redundant datasets can be stored in only a fraction of its origional memory
**/
class DemultiplexedReads
{
    public:

        // Constructor that initializes both vectors with x elements
        DemultiplexedReads(size_t readNum) 
        {
            lineBatch.reserve(readNum);
        }

        void store_demultiplexed_read(DemultiplexedLine& line)
        {
            lineBatch.push_back(line);;
        }

        size_t size() const
        {
            return lineBatch.size();
        }

        const DemultiplexedLine at(const int& i)
        {
            return(lineBatch.at(i));
        }
        const std::vector<DemultiplexedLine>& get_all_reads() const
        {
            return(lineBatch);
        }
        bool contains_DNA()
        {
            return containsDNA;
        }

    private:
        std::vector<DemultiplexedLine> lineBatch;
        bool containsDNA = false;
};

typedef std::shared_ptr<DemultiplexedReads> DemultiplexedReadsPtr; 
