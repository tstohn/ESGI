#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "dataTypes.hpp"

/**
 * @brief A line in the demultiplexed data.
 * Storing the UMI sequence, Antibody sequence, 
 * unique cell sequences (basically all barcode sequences added to each other) 
 * and a sequence for the treatment
 */
struct dataLine
{
    const char* umiSeq;
    const char* abName;
    const char* scID;
    const char* treatmentName;
};
typedef std::shared_ptr<dataLine> dataLinePtr;

/**
 * @brief A class storing all the demultiplexed barcodes.
 */
class UnprocessedDemultiplexedData
{
    public:

        UnprocessedDemultiplexedData()
        {
            uniqueChars = std::make_shared<UniqueCharSet>();
        }

        // add a dataLines to the vector
        void add(std::string& umiStr, std::string& abStr, std::string& singleCellStr, std::string& treatment)
        {
            //get unique pointer for all three strings
            dataLine line;
            line.umiSeq = uniqueChars->getUniqueChar(umiStr.c_str());
            line.abName = uniqueChars->getUniqueChar(abStr.c_str());
            line.scID = uniqueChars->getUniqueChar(singleCellStr.c_str());
            line.treatmentName = uniqueChars->getUniqueChar(treatment.c_str());

            //make a dataLinePtr from those unique strings
            dataLinePtr linePtr(std::make_shared<dataLine>(line));

            //add it to our dataStructure (3 entries have to be set)
            addDataLine(linePtr);
        }

        inline std::shared_ptr<UniqueCharSet> getUniqueBarcodes() const
        {
            return uniqueChars;
        }
        //return functions for our data, based on positions, UMI or AB/ SC barcodes
        /*inline const std::vector<dataLinePtr> getData() const
        {
            return data;
        }
        inline dataLinePtr getDataAt(int pos) const
        {
            return data.at(pos);
        }*/
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
        inline std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> getUniqueSc() const
        {
            return positionsOfSingleCell;
        }
        //we have to move the dataLine from the vector of lines for of oldUmi to newUmi
        //additionally inside this line we must update the new UMI sequence
        inline void changeUmi(const char* oldUmi, const char* newUmi, dataLinePtr oldLine)
        {
            //move old line from 'wrong' UMI key to right key
            remove(positionsOfUmi.at(oldUmi).begin(), positionsOfUmi.at(oldUmi).end(), oldLine);
            positionsOfUmi.at(newUmi).push_back(oldLine);
            oldLine->umiSeq = newUmi;
        }
        inline void setTreatmentDict(std::unordered_map<std::string, std::string > dict)
        {
            treatmentDict = dict;
        }
        inline void setProteinDict(std::unordered_map<std::string, std::string > dict)
        {
            proteinDict = dict;
        }
        inline std::string getProteinName(std::string barcode)
        {
            return proteinDict[barcode];
        }
        inline std::string getTreatmentName(std::string barcode)
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
            //we do not need to keep all lines in order, we r only interested in lines
            //for a UMI, AbSc, or a Sc
            //data.push_back(line);

            //2.) INSERT UMI POSITIONS
            //if umi already exists also store this new position
            if(positionsOfUmi.find(line->umiSeq) == positionsOfUmi.end())
            {
                std::vector<dataLinePtr> vec;
                vec.push_back(line);
                positionsOfUmi.insert(std::make_pair(line->umiSeq, vec));
            }
            //if not add this new umi with this actual position to map
            else
            {
                positionsOfUmi[line->umiSeq].push_back(line);
            }

            // 3.) INSERT ABSC POSITIONS
            //same for AbSingleCell
            std::string abScIdxStr = std::string((line->abName)) + std::string((line->scID));
            
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

            // 3.) INSERT SC POSITIONS
            std::string singleCell = std::string((line->scID));
            const char* singleCellChar = uniqueChars->getUniqueChar(singleCell.c_str());
            if(positionsOfSingleCell.find(singleCellChar) == positionsOfSingleCell.end())
            {
                std::vector<dataLinePtr> vec;
                vec.push_back(line);
                positionsOfSingleCell.insert(std::make_pair(singleCellChar, vec));
            }
            else
            {
                positionsOfSingleCell[singleCellChar].push_back(line);
            }
        }

        std::string uniqueAbSingleCellId(std::string abId, std::string singleCellId)
        {
            std::string id = abId + singleCellId;
            return id;
        }

        //hash tables storing all the positions of dataLines for a unique
        // a) UMI  b) SingleCell-AB combination
        std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> positionsOfUmi;
        std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> positonsOfABSingleCell;
        std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> positionsOfSingleCell;

        //all the string inside this class are stored only once, 
        //for all strings scID, Ab-name, treatment-name we store the string only once, and then ptrs to it
        std::shared_ptr<UniqueCharSet> uniqueChars;

        //dictionaries to map a barcode-sequence to the treatment, and Protein, 
        //those dicts are used in the very beginning when lines r parsed, so the real barcode sequence is never stored
        std::unordered_map<std::string, std::string > treatmentDict;
        std::unordered_map<std::string, std::string > proteinDict;
};