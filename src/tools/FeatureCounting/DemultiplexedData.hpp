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
    //variables that r set directly
    const char* umiSeq;
    const char* abName;
    const char* scID;
    const char* treatmentName;

    //variables set later on
    //class name and umi count are set when removing non-unique UMI read and collapsing the umis
    const char* cellClassname;
    unsigned long long umiCount = 0;
};
typedef std::shared_ptr<dataLine> umiDataLinePtr;
typedef std::shared_ptr<const dataLine> dataLinePtr;

//less operator to compare two dataLines, compares the length distance of the lines UMIs to the
//origional length that this line is supposed to have
struct less_than_umi
{
    less_than_umi(const int& origionalLength, 
    std::shared_ptr<std::unordered_map<const char*, std::vector<umiDataLinePtr>, 
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

        unsigned long long readNum1 = umiToReadsMap->at(line1->umiSeq).size();
        unsigned long long readNum2 = umiToReadsMap->at(line2->umiSeq).size();

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
    std::shared_ptr<std::unordered_map<const char*, std::vector<umiDataLinePtr>, 
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
            positionsOfUmiPtr = std::make_shared< std::unordered_map<const char*, std::vector<umiDataLinePtr>, CharHash, CharPtrComparator> >();
            positonsOfABSingleCellPtr = std::make_shared< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> >();
        }

        //during parsing of the demultiplexed reads we add all reads firstly to a tmp structure if we also parse guide reads
        //this functions stores the AB count data in a temporary vector
        void add_to_tmp_dataLines(std::string& umiStr, std::string& abStr, std::string& singleCellStr, std::string& treatment,
        std::vector<dataLine> abDataLines)
        {
            dataLine line;
            line.umiSeq = uniqueChars->getUniqueChar(umiStr.c_str());
            line.abName = uniqueChars->getUniqueChar(abStr.c_str());
            line.scID = uniqueChars->getUniqueChar(singleCellStr.c_str());
            line.treatmentName = uniqueChars->getUniqueChar(treatment.c_str());

            abDataLines.push_back(line);
        }
        //this fucntion stores the guide reads in a map, mapping scIds to the occurence of the different class labels
        void add_tmp_class_line(std::string& className, std::string& scId,
                    std::unordered_map< const char*, std::unordered_map< const char*, UnorderedSetCharPtr>>& scClasseCountDict,
                    const char* tmpUmi)
        {
            const char* uniqueUmi = uniqueChars->getUniqueChar(tmpUmi); //adding the temporary char* of the parsed line into our unique char dict

            const char* scCharPtr = uniqueChars->getUniqueChar(scId.c_str());
            const char* nameCharPtr = uniqueChars->getUniqueChar(className.c_str());

            //increase the count of this class for the specific cell
            if(scClasseCountDict.find(scCharPtr) == scClasseCountDict.end())
            {
                std::unordered_map<const char*, UnorderedSetCharPtr> map;
                UnorderedSetCharPtr set;
                set.insert(uniqueUmi);
                map.insert(std::make_pair(nameCharPtr, set));
                scClasseCountDict.insert(std::make_pair(scCharPtr, map));
            }
            //if not add this new umi with this actual position to map
            else
            {
                std::unordered_map<const char*, UnorderedSetCharPtr> innerMap = scClasseCountDict.at(scCharPtr);
                if(innerMap.find(nameCharPtr) == innerMap.end())
                {
                    std::unordered_map<const char*, UnorderedSetCharPtr> map;
                    UnorderedSetCharPtr set;
                    set.insert(uniqueUmi);
                    scClasseCountDict.at(scCharPtr).insert(std::make_pair(nameCharPtr, set));
                }
                else
                {
                    (scClasseCountDict.at(scCharPtr).at(nameCharPtr)).insert(uniqueUmi);
                }
            }
        }

        // add a dataLines to the vector
        void add_to_umiDict(const char* umiChar, std::string& abStr, std::string& singleCellStr, std::string& treatment)
        {
            //get unique pointer for all three strings
            dataLine line;
            line.umiSeq = uniqueChars->getUniqueChar(umiChar);
            line.abName = uniqueChars->getUniqueChar(abStr.c_str());
            line.scID = uniqueChars->getUniqueChar(singleCellStr.c_str());
            line.treatmentName = uniqueChars->getUniqueChar(treatment.c_str());

            //make a dataLinePtr from those unique strings
            umiDataLinePtr linePtr(std::make_shared<dataLine>(line));
            //add it to our dataStructure (3 entries have to be set)
            add_dataLine_to_umiDict(linePtr);
        }

        // add a dataLines to the vector
        void add_to_scAbDict(const char* umiChar, std::string& abStr, std::string& singleCellStr, std::string& treatment)
        {
            //get unique pointer for all three strings
            dataLine line;
            line.umiSeq = uniqueChars->getUniqueChar(umiChar);
            line.abName = uniqueChars->getUniqueChar(abStr.c_str());
            line.scID = uniqueChars->getUniqueChar(singleCellStr.c_str());
            line.treatmentName = uniqueChars->getUniqueChar(treatment.c_str());

            //make a dataLinePtr from those unique strings
            dataLinePtr linePtr(std::make_shared<dataLine>(line));

            //add it to our dataStructure (3 entries have to be set)
            add_dataLine_to_scabDict(linePtr);
        }

        // add a umi dataline to its final ABSc Dict structure (add a class name and add to dict)
        void add_to_scAbDict(const umiDataLinePtr& line, unsigned long long& umiCount, const char* cellClass = nullptr)
        {
            //add class name
            line->cellClassname = cellClass;
            line->umiCount = umiCount;
            
            //add it to our dataStructure (3 entries have to be set)
            add_dataLine_to_scabDict(line);
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
        inline std::vector<umiDataLinePtr> getDataWithUmi(const char* umi) const
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
        inline const std::shared_ptr< std::unordered_map<const char*, std::vector<umiDataLinePtr>, CharHash, CharPtrComparator>>getUniqueUmis() const
        {
            return positionsOfUmiPtr;
        }
        inline const std::shared_ptr< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator>> getUniqueAbSc() const
        {
            return positonsOfABSingleCellPtr;
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

            //std::string singleCell = std::string((oldLine->scID));
            //const char* singleCellChar = uniqueChars->getUniqueChar(singleCell.c_str());
            
            //remove all the old unnecessary lines from the dicts
            remove(positionsOfUmiPtr->at(oldUmi).begin(), positionsOfUmiPtr->at(oldUmi).end(), oldLine);
            remove(positonsOfABSingleCellPtr->at(abScIdxChar).begin(), positonsOfABSingleCellPtr->at(abScIdxChar).end(), oldLine);

            //generate the new line
            dataLine newLine = *oldLine;
            newLine.umiSeq = newUmi;
            std::shared_ptr<const dataLine> newLinePtr = std::make_shared<const dataLine>(newLine);
            std::shared_ptr<dataLine> newLineUmiPtr = std::make_shared<dataLine> (newLine);

            //and add the new line to all dicts
            positionsOfUmiPtr->at(newUmi).push_back(newLineUmiPtr);
            positonsOfABSingleCellPtr->at(abScIdxChar).push_back(newLinePtr);
        }
        inline void setTreatmentDict(std::unordered_map<std::string, std::string > dict)
        {
            treatmentDict = dict;
        }
        inline void setProteinDict(std::unordered_map<std::string, std::string > dict)
        {
            proteinDict = dict;
            if(proteinDict.empty())
            {
                mapFeatureNames = false;
            }
        }
        inline void setClassDict(std::unordered_map<std::string, std::string > dict)
        {
            classDict = dict;
        }
        inline void set_cell_to_class_dict(std::unordered_map< const char*, const char*> tmpCcClassMap)
        {
            scClassMap = tmpCcClassMap;
        }

        inline std::string getFeatureName(const std::string& barcode) const
        {
            if(mapFeatureNames == false){return barcode;}
            //if we have to map names
            if ( proteinDict.find(barcode) == proteinDict.end() ) 
            {
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
        inline const char* get_sc_class_name(const char* sc)
        {
            if(scClassMap.count(sc))
            {
                return(scClassMap.at(sc));
            }
            return(nullptr);
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

    private:

        void add_dataLine_to_scabDict(const dataLinePtr& line)
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
                positonsOfABSingleCellPtr->at(abScIdxChar).push_back(line);
            }
        }

        void add_dataLine_to_umiDict(const umiDataLinePtr& line)
        {
            if(positionsOfUmiPtr->find(line->umiSeq) == positionsOfUmiPtr->end())
            {
                std::vector<umiDataLinePtr> vec;
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
        //the first dict to store all demultiplexed lines (they r ordered by the umi, bcs first we do a umi sanity check - one umi == one sc/ab - and collapse umis)
        std::shared_ptr< std::unordered_map<const char*, std::vector<umiDataLinePtr>, CharHash, CharPtrComparator> > positionsOfUmiPtr;
        //after collapsing the umis we store all reads as const dataLines, we only still perform a umiMM corrections step within reads of each AB/SC 
        std::shared_ptr< std::unordered_map<const char*, std::vector<dataLinePtr>, CharHash, CharPtrComparator> > positonsOfABSingleCellPtr;

        std::unordered_map< const char*, const char*> scClassMap;

        //all the string inside this class are stored only once, 
        //for all strings scID, Ab-name, treatment-name we store the string only once, and then ptrs to it
        std::shared_ptr<UniqueCharSet> uniqueChars;

        //dictionaries to map a barcode-sequence to the treatment, and Protein, class
        //those dicts are used in the very beginning when lines r parsed, so the real barcode sequence is never stored
        std::unordered_map<std::string, std::string > treatmentDict;
        bool mapFeatureNames = true; //if FeatureCounting shall map feature barcodes to names or not
        //e.g.: for scRNAseq we assign gene names previously and do not need to map barcodes to names, for AB-barcodes this might still be necessary
        std::unordered_map<std::string, std::string > proteinDict;
        std::unordered_map<std::string, std::string > classDict;
};