#include <iostream>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include "BarcodeProcessingHandler.hpp"
#include "helper.hpp"

struct VectorHasher {
    int operator()(const std::vector<int> &V) const {
        int hash = V.size();
        for(auto &i : V) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

class umiQualityStat
{
    public:
        // adds data to the errorMap
        //input is following map: barcodeID(AB,BC1, bC2, ...) => vector of different barcodes found
        void add_value(const std::unordered_map<std::string, std::vector<std::string> >& barcodesForUmi);
        std::vector<std::string> get_keys_in_order()
        {
            return(barcodeOrder);
        }
        std::unordered_map< std::vector<int>, unsigned long long, VectorHasher> get_error_map()
        {
            return(errorMap);
        }

    private:
        //Map to store where we have several barcodes for a UMI
        // BC1 | BC2 | ABname | treatmentname => count
        // 1   | 1   | 1      | 1             => 92236
        // 1   | 2   | 1      | 1             => 100
        //e.g. we have many UMIs (92236) where the Sc-AB-treatment is truly unique
        //and 100 UMIs that have two different barcodes in BC-round 2
        std::unordered_map< std::vector<int>, unsigned long long, VectorHasher> errorMap;
        std::vector<std::string> barcodeOrder; // order of barcode types (AB, BC1, BC2, ...) in which the counts are stored in vector<int>
        //above in the errorMap
        std::mutex errorMapUpdateLock;
};

class UmiQuality
{
    public:

        UmiQuality(const BarcodeProcessingHandler& handler)
        {
            rawData = handler.getRawData();
        }
        //run the quality check: calls 1.) checkUniquenessOfUmis 2.) writeUmiQualityData
        void runUmiQualityCheck(const int& thread, const std::string& output);

    private:
    //private functions called in runUmiQualityCheck
        void checkUniquenessOfUmis(const std::vector<umiDataLinePtr>& uniqueUmis);
        void writeUmiQualityData(std::string output);

        //Statistic about the UMI quality
        umiQualityStat umiQualStat;
        //raw data: storing all demultiplexed dataLines
        UnprocessedDemultiplexedData rawData;
};