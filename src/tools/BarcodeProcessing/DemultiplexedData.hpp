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

    const char* cellClassname;
};
typedef std::shared_ptr<const dataLine> dataLinePtr;

//less operator to compare thwo dataLines, compares the length distance of the lines UMIs to the
//origional length that this line is supposed to have
struct less_than_umi
{
    less_than_umi(const int& origionalLength, 
    std::shared_ptr<std::unordered_map<const char*, std::vector<dataLinePtr>, 
                    CharHash, CharPtrComparator>>& umiToReadsMap)
    {
        this->origionalLength = origionalLength;
        this->umiToReadsMap = umiToReadsMap;
    }
    inline bool operator() (const dataLinePtr& line1, const dataLinePtr& line2)
    {
        int lenDiff1 = std::strlen(line1->umiSeq) - origionalLength;
        lenDiff1 = sqrt(lenDiff1*lenDiff1);
        int lenDiff2 = std::strlen(line2->umiSeq) - origionalLength;
        lenDiff2 = sqrt(lenDiff2*lenDiff2);

        unsigned long long readNum1 = (*umiToReadsMap)[line1->umiSeq].size();
        unsigned long long readNum2 = (*umiToReadsMap)[line2->umiSeq].size();

        if(lenDiff1 != lenDiff2)
        {
            return (lenDiff1 < lenDiff2);
        }
        else
        {
            return(readNum1 > readNum2);
        }
    }
    int origionalLength;
    std::shared_ptr<std::unordered_map<const char*, std::vector<dataLinePtr>, 
                    CharHash, CharPtrComparator>> umiToReadsMap;
};

/**
 * @brief A class storing all the demultiplexed barcodes.
 */
class UnprocessedDemultiplexedData
{
    public:

        UnprocessedDemultiplexedData()
        {
            uniqueChars = std::make_shared<UniqueCharSet>();
            positionsOfUmiPtr = std::make_shared< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> >();
            positonsOfABSingleCellPtr = std::make_shared< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> >();
            positionsOfSingleCellPtr = std::make_shared< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> >();
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

                // add a dataLines to the vector
        void add_to_umiDict(std::string& umiStr, std::string& abStr, std::string& singleCellStr, std::string& treatment)
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
            add_dataLine_to_umiDict(linePtr);
        }

                        // add a dataLines to the vector
        void add_to_scAbDict(std::string& umiStr, std::string& abStr, std::string& singleCellStr, std::string& treatment)
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
            add_dataLine_to_scabDict(linePtr);
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
            return positionsOfUmiPtr->at(umi);
        }
        inline std::vector<dataLinePtr> getDataWithAbSc(const char* ab, const char* sc) const
        {
            std::string abScIdxStr = std::string((ab)) + std::string((sc));
            return positonsOfABSingleCellPtr->at(abScIdxStr.c_str());
        }
        inline std::vector<dataLinePtr> getDataWithAbSc(const char* absc) const
        {
            return positonsOfABSingleCellPtr->at(absc);
        }
        inline const std::shared_ptr< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator>>getUniqueUmis() const
        {
            return positionsOfUmiPtr;
        }
        inline const std::shared_ptr< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator>> getUniqueAbSc() const
        {
            return positonsOfABSingleCellPtr;
        }
        inline const std::shared_ptr< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator>> getUniqueSc() const
        {
            return positionsOfSingleCellPtr;
        }
        //we have to move the dataLine from the vector of lines for of oldUmi to newUmi
        //additionally inside this line we must update the new UMI sequence
        //THIS FUNCTION IS NOT TESTED, also not used at the moment: be aware when using !!!!
        inline void changeUmi(const char* oldUmi, const char* newUmi, dataLinePtr oldLine)
        {
            //we have 3 dictionaries to update

            //we have the key for the umiDict, but need the two keys for the other dicts
            std::string abScIdxStr = std::string((oldLine->abName)) + std::string((oldLine->scID));
            const char* abScIdxChar = uniqueChars->getUniqueChar(abScIdxStr.c_str());

            std::string singleCell = std::string((oldLine->scID));
            const char* singleCellChar = uniqueChars->getUniqueChar(singleCell.c_str());
            
            //remove all the old unnecessary lines from the dicts
            remove(positionsOfUmiPtr->at(oldUmi).begin(), positionsOfUmiPtr->at(oldUmi).end(), oldLine);
            remove(positonsOfABSingleCellPtr->at(abScIdxChar).begin(), positonsOfABSingleCellPtr->at(abScIdxChar).end(), oldLine);
            remove(positionsOfSingleCellPtr->at(singleCellChar).begin(), positionsOfSingleCellPtr->at(singleCellChar).end(), oldLine);

            //generate the new line
            dataLine newLine = *oldLine;
            newLine.umiSeq = newUmi;
            std::shared_ptr<const dataLine> newLinePtr = std::make_shared<const dataLine>(newLine);
            
            //and add the new line to all dicts
            positionsOfUmiPtr->at(newUmi).push_back(newLinePtr);
            positonsOfABSingleCellPtr->at(abScIdxChar).push_back(newLinePtr);
            positionsOfSingleCellPtr->at(singleCellChar).push_back(newLinePtr);
        }
        inline void setTreatmentDict(std::unordered_map<std::string, std::string > dict)
        {
            treatmentDict = dict;
        }
        inline void setProteinDict(std::unordered_map<std::string, std::string > dict)
        {
            proteinDict = dict;
        }
        inline void setClassDict(std::unordered_map<std::string, std::string > dict)
        {
            classDict = dict;
        }
        inline std::string getProteinName(const std::string& barcode) const
        {
             if ( proteinDict.find(barcode) == proteinDict.end() ) {
                std::cerr << "Barcode is not in protein dict, check protein barcodes, protein ID and protein File:\n" << barcode << "\n";
                exit(EXIT_FAILURE);
            }
            return proteinDict.at(barcode);
        }
        inline std::string getTreatmentName(const std::string& barcode) const
        {
            if ( treatmentDict.find(barcode) == treatmentDict.end() ) {
                std::cerr << "Barcode is not in treatment dict, check treatment barcodes, treatmend ID and treatment File:\n" << barcode << "\n";
                exit(EXIT_FAILURE);
            }
            return treatmentDict.at(barcode);
        }
        inline std::string getClassName(const std::string& barcode) const
        {
            if ( classDict.find(barcode) == classDict.end() ) {
                std::cerr << "Barcode is not in Class dict, check Class barcodes, Class ID and Class File:\n" << barcode << "\n";
                exit(EXIT_FAILURE);
            }
            return classDict.at(barcode);
        }
        inline std::string get_protein_or_class_name(const std::string& barcode, bool& is_class) const
        {
            if ( proteinDict.find(barcode) == proteinDict.end() ) 
            {
                std::string className = getClassName(barcode);
                is_class = true;
                return className;
            }
            //else
            return proteinDict.at(barcode);
        }
        inline bool check_class() const
        {
            if(classDict.empty())
            {
                return false;
            }
            //else
            return true;
        }
        void add_class_line(std::string& className, std::string& scId)
        {
            const char* scCharPtr = uniqueChars->getUniqueChar(scId.c_str());
            const char* nameCharPtr = uniqueChars->getUniqueChar(className.c_str());

            //increase the count of this class for the specific cell
            if(scClasseCountDict.find(scCharPtr) == scClasseCountDict.end())
            {
                std::unordered_map<const char*, unsigned long long> map;
                map.insert(std::make_pair(nameCharPtr, 1));
                scClasseCountDict.insert(std::make_pair(scCharPtr, map));
            }
            //if not add this new umi with this actual position to map
            else
            {
                std::unordered_map<const char*, unsigned long long> innerMap = scClasseCountDict.at(scCharPtr);
                if(innerMap.find(nameCharPtr) == innerMap.end())
                {
                    scClasseCountDict.at(scCharPtr).insert(std::make_pair(nameCharPtr, 1));
                }
                else
                {
                    ++(scClasseCountDict.at(scCharPtr).at(nameCharPtr));
                }
            }
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
            if(positionsOfUmiPtr->find(line->umiSeq) == positionsOfUmiPtr->end())
            {
                std::vector<dataLinePtr> vec;
                vec.push_back(line);
                positionsOfUmiPtr->insert(std::make_pair(line->umiSeq, vec));
            }
            //if not add this new umi with this actual position to map
            else
            {
                (*positionsOfUmiPtr)[line->umiSeq].push_back(line);
            }

            // 3.) INSERT ABSC POSITIONS
            //same for AbSingleCell
            std::string abScIdxStr = std::string((line->abName)) + std::string((line->scID));
            
            const char* abScIdxChar = uniqueChars->getUniqueChar(abScIdxStr.c_str());
            if(positonsOfABSingleCellPtr->find(abScIdxChar) == positonsOfABSingleCellPtr->end())
            {
                std::vector<dataLinePtr> vec;
                vec.push_back(line);
                positonsOfABSingleCellPtr->insert(std::make_pair(abScIdxChar, vec));
            }
            //if not add this new umi with this actual position to map
            else
            {
                (*positonsOfABSingleCellPtr)[abScIdxChar].push_back(line);
            }

            // 4.) INSERT SC POSITIONS
            std::string singleCell = std::string((line->scID));
            const char* singleCellChar = uniqueChars->getUniqueChar(singleCell.c_str());
            if(positionsOfSingleCellPtr->find(singleCellChar) == positionsOfSingleCellPtr->end())
            {
                std::vector<dataLinePtr> vec;
                vec.push_back(line);
                positionsOfSingleCellPtr->insert(std::make_pair(singleCellChar, vec));
            }
            else
            {
                (*positionsOfSingleCellPtr)[singleCellChar].push_back(line);
            }
        }

        void add_dataLine_to_scabDict(dataLinePtr line)
        {
         // 3.) INSERT ABSC POSITIONS
            //same for AbSingleCell
            std::string abScIdxStr = std::string((line->abName)) + std::string((line->scID));
            
            const char* abScIdxChar = uniqueChars->getUniqueChar(abScIdxStr.c_str());
            if(positonsOfABSingleCellPtr->find(abScIdxChar) == positonsOfABSingleCellPtr->end())
            {
                std::vector<dataLinePtr> vec;
                vec.push_back(line);
                positonsOfABSingleCellPtr->insert(std::make_pair(abScIdxChar, vec));
            }
            //if not add this new umi with this actual position to map
            else
            {
                (*positonsOfABSingleCellPtr)[abScIdxChar].push_back(line);
            }

            // 4.) INSERT SC POSITIONS
            std::string singleCell = std::string((line->scID));
            const char* singleCellChar = uniqueChars->getUniqueChar(singleCell.c_str());
            if(positionsOfSingleCellPtr->find(singleCellChar) == positionsOfSingleCellPtr->end())
            {
                std::vector<dataLinePtr> vec;
                vec.push_back(line);
                positionsOfSingleCellPtr->insert(std::make_pair(singleCellChar, vec));
            }
            else
            {
                (*positionsOfSingleCellPtr)[singleCellChar].push_back(line);
            }
        }

        void add_dataLine_to_umiDict(dataLinePtr line)
        {
            if(positionsOfUmiPtr->find(line->umiSeq) == positionsOfUmiPtr->end())
            {
                std::vector<dataLinePtr> vec;
                vec.push_back(line);
                positionsOfUmiPtr->insert(std::make_pair(line->umiSeq, vec));
            }
            //if not add this new umi with this actual position to map
            else
            {
                (*positionsOfUmiPtr)[line->umiSeq].push_back(line);
            }
        }

        std::string uniqueAbSingleCellId(std::string abId, std::string singleCellId)
        {
            std::string id = abId + singleCellId;
            return id;
        }

        //hash tables storing all the positions of dataLines for a unique (stored ad shared_ptr so handing those maps over during processing is cheap)
        // a) UMI  b) SingleCell-AB combination
        std::shared_ptr< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> > positionsOfUmiPtr;
        std::shared_ptr< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> > positonsOfABSingleCellPtr;
        std::shared_ptr< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> > positionsOfSingleCellPtr;

        //all the string inside this class are stored only once, 
        //for all strings scID, Ab-name, treatment-name we store the string only once, and then ptrs to it
        std::shared_ptr<UniqueCharSet> uniqueChars;

        //dictionaries to map a barcode-sequence to the treatment, and Protein, class
        //those dicts are used in the very beginning when lines r parsed, so the real barcode sequence is never stored
        std::unordered_map<std::string, std::string > treatmentDict;
        std::unordered_map<std::string, std::string > proteinDict;
        std::unordered_map<std::string, std::string > classDict;
};