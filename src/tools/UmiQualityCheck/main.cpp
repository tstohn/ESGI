#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

#include "BarcodeProcessingHandler.hpp"
#include "UmiQualityHelper.hpp"

using namespace boost::program_options;

/**
 * Basic Tool to run a Quality check on the UMIs.
 * It generates an Overview in which barcode-pattern we see most
 * different barcodes for a unique UMI 
 * (useful to analyze random barcode combinations in CI experiments)
 * 
 * Input is similar to the BarcodeProcessing tool.
 * Output is a file of following format:
 *  Map to store where we have several barcodes for a UMI
 *    BC1 | BC2 | ABname | treatmentname => count
 *    1   | 1   | 1      | 1             => 92236
 *    1   | 2   | 1      | 1             => 100
 *   e.g. we have many UMIs (92236) where the Sc-AB-treatment is truly unique
 *   and 100 UMIs that have two different barcodes in BC-round 2
 **/

bool parse_arguments(char** argv, int argc, std::string& inFile,  std::string& outFile, int& threats, 
                     std::string& barcodeFile, std::string& barcodeIndices, int& umiMismatches,
                     std::string& abFile, int& abIdx, std::string& treatmentFile, int& treatmentIdx)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input,i", value<std::string>(&inFile)->required(), "directory of files or single file in fastq(.gz) format")
            ("output,o", value<std::string>(&outFile)->required(), "output file")

            ("barcodeList,b", value<std::string>(&(barcodeFile)), "file with a list of all allowed well barcodes (comma seperated barcodes across several rows)\
            the row refers to the correponding bracket enclosed sequence substring. E.g. for two bracket enclosed substrings in out sequence a possible list could be:\
            AGCTTCGAG,ACGTTCAGG\nACGTCTAGACT,ATCGGCATACG,ATCGCGATC,ATCGCGCATAC. This can be the same list as it was for FastqParser.")
            ("antibodyList,a", value<std::string>(&(abFile)), "file with a list of all antbodies used, should be in same order as the ab-barcodes in the barcodeList.")
            ("antibodyIndex,x", value<int>(&abIdx), "Index used for antibody distinction.")
            ("groupList,g", value<std::string>(&(treatmentFile)), "file with a list of all groups (e.g.treatments) used, should be in same order as the specific arcodes in the barcodeList. \
            If tis argument is given, you must also add the index of barcodes used for grouping")
            ("GroupingIndex,y", value<int>(&treatmentIdx), "Index used to group cells(e.g. by treatment). This is the x-th barcode from the barcodeFile (0 indexed).")

            ("CombinatorialIndexingBarcodeIndices,c", value<std::string>(&(barcodeIndices))->required(), "comma seperated list of indexes, that are used during \
            combinatorial indexing and should distinguish a unique cell. Be aware that this is the index of the line inside the barcodeList file (see above). \
            This file ONLY includes lines for the varying sequences (except UMI). Therefore the index is not the same as the position in the whole sequence \
            if constant or UMI-seq are present. Index starts with zero.")
            
            ("mismatches,u", value<int>(&umiMismatches)->default_value(2), "number of allowed mismatches in a UMI. The nucleotides in the beginning and end do NOT count.\
            Since the UMI is defined as the sequence between the last and first match of neighboring sequences, bases of mismatches could be in the beginning/ end.")
            ("thread,t", value<int>(&threats)->default_value(5), "number of threads")
            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";
            std::cout << "EXAMPLE CALL:\n ./bin/umiqual -i <inFile> ... \n";
            return false;
        }

        notify(vm);
    }
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    return true;
}

// generate a dictionary to map sequences to AB(proteins)
std::unordered_map<std::string, std::string > generateProteinDict(std::string abFile, int abIdx, 
                                                                                              const std::vector<std::string>& abBarcodes)
{
    std::unordered_map<std::string, std::string > map;
    std::vector<std::string> proteinNames;

    std::ifstream abFileStream(abFile);
    for(std::string line; std::getline(abFileStream, line);)
    {
        std::string delimiter = ",";
        std::string seq;
        size_t pos = 0;
        std::vector<std::string> seqVector;
        while ((pos = line.find(delimiter)) != std::string::npos) 
        {
            seq = line.substr(0, pos);
            line.erase(0, pos + 1);
            proteinNames.push_back(seq);
        }
        seq = line;
        proteinNames.push_back(seq);
    }

    assert(abBarcodes.size() == proteinNames.size());
    for(int i = 0; i < abBarcodes.size(); ++i)
    {
        map.insert(std::make_pair(abBarcodes.at(i), proteinNames.at(i)));
    }
    abFileStream.close();

    return map;
}

// generate a dictionary to map sequences to treatments
std::unordered_map<std::string, std::string > generateTreatmentDict(std::string treatmentFile, int treatmentIdx,
                                                                                                const std::vector<std::string>& treatmentBarcodes)
{
    std::unordered_map<std::string, std::string > map;
    std::vector<std::string> treatmentNames;

    std::ifstream treatmentFileStream(treatmentFile);
    for(std::string line; std::getline(treatmentFileStream, line);)
    {
        std::string delimiter = ",";
        std::string seq;
        size_t pos = 0;
        std::vector<std::string> seqVector;
        while ((pos = line.find(delimiter)) != std::string::npos) 
        {
            seq = line.substr(0, pos);
            line.erase(0, pos + 1);
            treatmentNames.push_back(seq);

        }
        seq = line;
        treatmentNames.push_back(seq);
    }
    assert(treatmentNames.size() == treatmentBarcodes.size());
    for(int i = 0; i < treatmentBarcodes.size(); ++i)
    {
        map.insert(std::make_pair(treatmentBarcodes.at(i), treatmentNames.at(i)));
    }
    treatmentFileStream.close();

    return map;
}


int main(int argc, char** argv)
{

    std::string inFile;
    std::string outFile;
    std::string barcodeFile;
    std::string barcodeIndices;
    int thread;
    int umiMismatches;

    //data for protein(ab) and treatment information
    std::string abFile; 
    int abIdx;
    std::string treatmentFile;
    int treatmentIdx;
    std::vector<std::string> abBarcodes;
    std::vector<std::string> treatmentBarcodes;

    parse_arguments(argv, argc, inFile, outFile, thread, barcodeFile, barcodeIndices, umiMismatches, abFile, abIdx, treatmentFile, treatmentIdx);
    
    //generate the dictionary of barcode alternatives to idx
    NBarcodeInformation barcodeIdData;
    generateBarcodeDicts(barcodeFile, barcodeIndices, barcodeIdData, abBarcodes, abIdx, &treatmentBarcodes, treatmentIdx);
    BarcodeProcessingHandler dataParser(barcodeIdData);
    //hack to prevent that the demultiplexed reads are written directly into the ScAb-matrix (which happens if the threshold to retains UMIs is set to 0)
    dataParser.setUmiFilterThreshold(-1.0);

    if(!abFile.empty())
    {
        std::unordered_map<std::string, std::string > map = generateProteinDict(abFile, abIdx, abBarcodes);
        dataParser.addProteinData(map);
    }
    if(!treatmentFile.empty())
    {
        std::unordered_map<std::string, std::string > map = generateTreatmentDict(treatmentFile, treatmentIdx, treatmentBarcodes);
        dataParser.addTreatmentData(map);
    }
    //parse the file of demultiplexed barcodes and
    //add all the data to the Unprocessed Demultiplexed Data (stored in rawData)
    // (AB, treatment already are mapped to their real names, scID is a concatenation of numbers for each barcode in
    //each abrcoding round, seperated by a dot)
    dataParser.parse_combined_file(inFile, thread);

    UmiQuality umiCheck(dataParser);
    umiCheck.runUmiQualityCheck(thread, outFile);

    return(EXIT_SUCCESS);
}