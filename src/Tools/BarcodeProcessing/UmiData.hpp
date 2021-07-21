#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <vector>
#include <thread>
#include <pthread.h>
#include <unordered_map>
#include <sstream>
#include <climits>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "dataTypes.hpp"

struct dataLine
{
    const char* umi_seq;
    const char* ab_seq;
    const char* cell_seq;
};
typedef std::shared_ptr<dataLine> dataLinePtr;
 
struct CIBarcode
{
    std::vector<std::unordered_map<std::string, int> > barcodeIdDict; // map of barcode to id, maps are in order of their occurence in the fastqRead
    // ids of the CI barcode within the barcode file (includes only barcodes for variable sequence regions)
    std::vector<int> ciBarcodeIndices; // more or less only of temporary usage, during generation of barcode map vector
};

class UmiData
{
    public:

        ~UmiData()
        {
            data.clear();
            positionsOfUmi.clear();
            posiitonsOfABSingleCell.clear();
            uniqueChars.~UniqueCharSet();
        }

        // add a dataLines to the vector
        void add(std::string& umiStr, std::string& abStr, std::string& singleCellStr)
        {
            //get unique pointer for all three string
            dataLine line;
            line.umi_seq = uniqueChars.getUniqueChar(umiStr.c_str());
            line.ab_seq = uniqueChars.getUniqueChar(abStr.c_str());
            line.cell_seq = uniqueChars.getUniqueChar(singleCellStr.c_str());

            //make a dataLinePtr from those unique string

            //add it to our dataStructure (3 entries have to be set)
            
        }

        //return all dataLines of this specific UMI
        std::vector<dataLinePtr> getDataLinesForUmi(const char* umi_seq)
        {
            std::vector<dataLinePtr> returnLines;

            return returnLines;
        }
        //return all dataLines of this specific AB-SingleCell combination
        std::vector<dataLinePtr> getDataLinesForABSingleCell(const char* abBarcode, const char* singleCellBarcode)
        {
            std::vector<dataLinePtr> returnLines;

            return returnLines;
        }

    private:

        std::vector<dataLinePtr> data;
        std::unordered_map<std::string, std::vector<int> > positionsOfUmi;
        std::unordered_map<std::string, std::vector<int> > posiitonsOfABSingleCell;

        //all the string inside this class are stored only once, 
        //set of all the unique barcodes we use, and we only pass pointers to those
        UniqueCharSet uniqueChars;

};

class UmiDataParser
{


    public:

        UmiDataParser(CIBarcode barcodeIdData) : barcodeDict(barcodeIdData){}
        ~UmiDataParser()
        {
            data.~UmiData();
        }

        void parseFile(const std::string fileName, const int& thread);

    private:

        //parse the file, store each line in our data structure
        void addFastqReadToUmiData(const std::string& line);
        void parseBarcodeLines(std::istream* instream, const int& totalReads, int& currentReads);

        //fill the CIBarcode dict: mapping of barcode alternatives to an unique idx
        void getCiBarcodeInWholeSequence(const std::string& line);
        //map all the barcodes of CI to a unique 'number' string as SingleCellIdx
        std::string generateSingleCellIndexFromBarcodes(std::vector<std::string> ciBarcodes);

        //data structure storing lines with: UMI, AB_id, SingleCell_id
        UmiData data;
        //Dictionary used to generate the dataLines, maps for each barcode in teh sequence all
        //possibilities to an idx
        CIBarcode barcodeDict;
        //
        std::vector<int> fastqReadBarcodeIdx; // ids of the CI barcode within the whole string of all barcodes
        int abIdx = INT_MAX;
        int umiIdx = INT_MAX;
};
