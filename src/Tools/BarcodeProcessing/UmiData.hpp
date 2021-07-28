#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "dataTypes.hpp"


struct abLine
{
    std::shared_ptr<std::string> ab_seq;
    std::shared_ptr<std::string> treatment;
    
    const char* cell_seq;
    int ab_cout = 0;
}; 

struct dataLine
{
    const char* umi_seq;
    const char* ab_seq;
    const char* cell_seq;
    const char* treatment_seq;
};
typedef std::shared_ptr<dataLine> dataLinePtr;

class UmiData
{
    public:

        UmiData()
        {
            uniqueChars = std::make_shared<UniqueCharSet>();
        }

        // add a dataLines to the vector
        void add(std::string& umiStr, std::string& abStr, std::string& singleCellStr, std::string& treatment)
        {
            //get unique pointer for all three string
            dataLine line;
            line.umi_seq = uniqueChars->getUniqueChar(umiStr.c_str());
            line.ab_seq = uniqueChars->getUniqueChar(abStr.c_str());
            line.cell_seq = uniqueChars->getUniqueChar(singleCellStr.c_str());
            line.treatment_seq = uniqueChars->getUniqueChar(treatment.c_str());

            //make a dataLinePtr from those unique string
            dataLinePtr linePtr(std::make_shared<dataLine>(line));

            //add it to our dataStructure (3 entries have to be set)

            addDataLine(linePtr);
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
        std::shared_ptr<UniqueCharSet> getUniqueBarcodes()
        {
            return uniqueChars;
        }

        //return functions for our data, based on positions, UMI or AB/ SC barcodes
        inline const std::vector<dataLinePtr> getData() const
        {
            return data;
        }
        inline dataLinePtr getDataAt(int pos) const
        {
            return data.at(pos);
        }
        inline std::vector<dataLinePtr> getDataWithUmi(const char* umi) const
        {
            return positionsOfUmi.at(umi);
        }
        inline std::vector<dataLinePtr> getDataWithAbSc(const char* ab, const char* sc) const
        {
            std::string abScIdxStr = std::string((ab)) + std::string((sc));
            return positonsOfABSingleCell.at(abScIdxStr.c_str());
        }
        inline std::vector<dataLinePtr> getDataWithAbSc(const char* absc) const
        {
            return positonsOfABSingleCell.at(absc);
        }
        inline std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> getUniqueUmis() const
        {
            return positionsOfUmi;
        }
        inline std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> getUniqueAbSc() const
        {
            return positonsOfABSingleCell;
        }
        inline void changeUmi(const char* oldUmi, const char* newUmi, dataLinePtr line)
        {
            remove(positionsOfUmi.at(oldUmi).begin(), positionsOfUmi.at(oldUmi).end(), line);
            positionsOfUmi.at(newUmi).push_back(line);
        }
        inline void setTreatmentDict(std::unordered_map<std::string, std::shared_ptr<std::string>> dict)
        {
            treatmentDict = dict;
        }
        inline void setProteinDict(std::unordered_map<std::string, std::shared_ptr<std::string>> dict)
        {
            proteinDict = dict;
        }
        inline std::shared_ptr<std::string> getProteinName(std::string barcode)
        {
            return proteinDict[barcode];
        }
        inline std::shared_ptr<std::string> getTreatmentName(std::string barcode)
        {
            if ( treatmentDict.find(barcode) == treatmentDict.end() ) {
                std::cerr << "Barcode is not in treatment dict, check treatment barcodes, treatmend ID and treatment File:\n" << barcode << "\n";
                exit(EXIT_FAILURE);
            }
            return treatmentDict[barcode];
        }

    private:

        void addDataLine(dataLinePtr line)
        {
            //1.) ADD LINE
            data.push_back(line);

            //2.) INSERT UMI POSITIONS
            //if umi already exists also store this new position
            if(positionsOfUmi.find(line->umi_seq) == positionsOfUmi.end())
            {
                std::vector<dataLinePtr> vec;
                vec.push_back(line);
                positionsOfUmi.insert(std::make_pair(line->umi_seq, vec));
            }
            //if not add this new umi with this actual position to map
            else
            {
                positionsOfUmi[line->umi_seq].push_back(line);
            }

            // 3.) INSERT ABSC POSITIONS
            //same for AbSingleCell
            std::string abScIdxStr = std::string((line->ab_seq)) + std::string((line->cell_seq));
            
            const char* abScIdxChar = uniqueChars->getUniqueChar(abScIdxStr.c_str());
            if(positonsOfABSingleCell.find(abScIdxChar) == positonsOfABSingleCell.end())
            {
                std::vector<dataLinePtr> vec;
                vec.push_back(line);
                positonsOfABSingleCell.insert(std::make_pair(abScIdxChar, vec));
            }
            //if not add this new umi with this actual position to map
            else
            {
                positonsOfABSingleCell[abScIdxChar].push_back(line);
            }
        }

        std::string uniqueAbSingleCellId(std::string abId, std::string singleCellId)
        {
            std::string id = abId + singleCellId;
            return id;
        }

        std::vector<dataLinePtr> data;
        std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> positionsOfUmi;
        std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> positonsOfABSingleCell;

        //all the string inside this class are stored only once, 
        //set of all the unique barcodes we use, and we only pass pointers to those
        std::shared_ptr<UniqueCharSet> uniqueChars;

        //dictionaries to map a barcode-sequence to the treatment, and Protein
        std::unordered_map<std::string, std::shared_ptr<std::string> > treatmentDict;
        std::unordered_map<std::string, std::shared_ptr<std::string> > proteinDict;
};