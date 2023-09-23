#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/version.hpp>

#include "BarcodeProcessingHandler.hpp"
using namespace boost::program_options;

/**
 * Tool to count antobodies for single cell in demultiplexed data.
 * It can consist of reads from ABs as well as guides. If guide reads r present two files with
 * the corresponding barcode sequences and their names (e.g. cell types) should be present. In this
 * case the barcodeList should still ONLY contain the AB barcodes.
 * 
 * If additional guide reads r present, the AB and guide reads can be in ONE file, or they can be in two seperate files, one for the AB
 * and one for the guide reads.
 * 
 * Output is a file with AB counts for all found single cells.
 * The UmiProcessing file contains the counts for each UMI that is in the end considered for the AB counts.
 * This is so to speak the number of PCR duplicates that we have for each single count of an AB within a single cell
 * This means guide reads are neglected, as well are all reads neglected that are filtered out during the processing (e.g. bcs of non unique UMI in reads, etc.)
 * 
 * Algorithm does following:
 * 1.) If guides are present, map each single cell to their guide, the guide must be unique for a single cell for >= 90% of the guide reads
 * 2.) compare all reads form a single UMI, and remove reads that are not unique (same AB/ single cell barcodes) and represent 90% of the reads
 * 3.) collapse all UMIs and also allow or mismatches between UMIs (only checked for reads of same AB and Single cell)
 * */

bool parse_arguments(char** argv, int argc, std::string& inFile,  std::string& outFile, int& threats, 
                     std::string& barcodeFile, std::string& barcodeIndices, int& umiMismatches,
                     std::string& abFile, int& abIdx, std::string& treatmentFile, int& treatmentIdx,
                     std::string& classSeqFile, std::string& classNameFile, double& umiThreshold,
                     bool& scClassConstraint, std::string& guideReadsFile, bool& umiRemoval,
                     std::string& umiSingleCellIdx)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input,i", value<std::string>(&inFile)->required(), "input file of demultiplexed reads for ABs in Single cells in tsv.gz format (input must be gzipped)")
            ("output,o", value<std::string>(&outFile)->required(), "output file with all split barcodes")

            ("barcodeList,b", value<std::string>(&(barcodeFile)), "file with a list of all allowed well barcodes (comma seperated barcodes across several rows)\
            the row refers to the correponding bracket enclosed sequence substring. E.g. for two bracket enclosed substrings in out sequence a possible list could be:\
            AGCTTCGAG,ACGTTCAGG\nACGTCTAGACT,ATCGGCATACG,ATCGCGATC,ATCGCGCATAC. This can be the same list as it was for FastqParser. Do not include the barcodes for the\
            guide reads here; add them as two seperate files guide reads are present.")
            ("antibodyList,a", value<std::string>(&(abFile)), "file with a list of all antbodies used, should be in same order as the ab-barcodes in the barcodeList.")
            ("antibodyIndex,x", value<int>(&abIdx), "Index used for antibody distinction. This is the x-th barcode from the barcodeFile (0 indexed)")
            ("groupList,d", value<std::string>(&(treatmentFile)), "file with a list of all groups (e.g.treatments) used, should be in same order as the specific arcodes in the barcodeList. \
            If this argument is given, you must also add the index of barcodes used for grouping")
            ("GroupingIndex,y", value<int>(&treatmentIdx), "Index used to group cells(e.g. by treatment). This is the x-th barcode from the barcodeFile (0 indexed).")

            ("classSeq,g", value<std::string>(&(classSeqFile)), "file with the sequences that define origin of cells (e.g. sgRNA sequences along the experiment)")
            ("className,n", value<std::string>(&(classNameFile)), "file with names to replace the sequence of origin")
            ("guideFile,j", value<std::string>(&guideReadsFile)->default_value(""), "file with all the demultiplexed guide reads.")

            ("CombinatorialIndexingBarcodeIndices,c", value<std::string>(&(barcodeIndices))->default_value(""), "comma seperated list of indexes, that are used during \
            combinatorial indexing and should distinguish a unique cell. Be aware that this is the index of the line inside the barcodeList file (see above). \
            This file ONLY includes lines for the varying sequences (except UMI). Therefore the index is not the same as the position in the whole sequence \
            if constant or UMI-seq are present. Index starts with zero.")

            ("mismatches,u", value<int>(&umiMismatches)->default_value(1), "number of allowed mismatches in a UMI. The nucleotides in the beginning and end do NOT count.\
            Since the UMI is defined as the sequence between the last and first match of neighboring sequences, bases of mismatches could be in the beginning/ end.")
            ("thread,t", value<int>(&threats)->default_value(5), "number of threads")
            ("umiThreshold,f", value<double>(&umiThreshold)->default_value(0.0), "threshold for filtering UMIs. E.g. if set to 0.9 we only retain reads of a UMI, if more \
            than 90percent of them have the same SC-AB combination. All other reads are deleted. Keep at 0 if UMIs should not be removed.")
            ("umiRemoval,z", value<bool>(&umiRemoval)->default_value(true), "Set to false if UMIs should NOT be collapsed.")
            ("scClassConstraint,k", value<bool>(&scClassConstraint)->default_value(true), "Boolean to store whether sc reads should be removed if we find no guide read for them. \
            If set to false reads for no guide are given the class wildtype.")
            ("umiSingleCellIdx,s", value<std::string>(&umiSingleCellIdx)->default_value(""), "In case the Single-Cell barcode(s) are not given beforehand.\
            In this case this Idx is the (0 indexed) position of the X-barcode, which should be used as single cell identifier.\
            E.g. for 10X when we have single cell indices that can not be distiguished from a UMI when mapping. This position must be a \
            sequence with the [X] pattern.")

            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";
            std::cout << "EXAMPLE CALL:\n ./bin/processing -i <inFile> ... \n";
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

// generate a dictionary to map sequences to treatments
std::unordered_map<std::string, std::string > generateClassDict(const std::string& classSeqFile,
                                                                const std::string& classNameFile)
{
    std::unordered_map<std::string, std::string > map;
    std::vector<std::string> names;
    std::vector<std::string> seqs;

    //parse all sequences
    std::ifstream seqFileStream(classSeqFile);
    for(std::string line; std::getline(seqFileStream, line);)
    {
        std::string delimiter = ",";
        std::string seq;
        size_t pos = 0;
        std::vector<std::string> seqVector;
        while ((pos = line.find(delimiter)) != std::string::npos) 
        {
            seq = line.substr(0, pos);
            line.erase(0, pos + 1);
            seqs.push_back(seq);
        }
        seq = line;
        seqs.push_back(seq);
    }
    seqFileStream.close();
    if(seqs.empty())
    {
        std::cout << "ERROR: Could not parse any sequence for guides! Check the guide sequence file.\n";
        exit(EXIT_FAILURE);
    }

    //parse all names
    std::ifstream nameFileStream(classNameFile);
    for(std::string line; std::getline(nameFileStream, line);)
    {
        std::string delimiter = ",";
        std::string seq;
        size_t pos = 0;
        std::vector<std::string> seqVector;
        while ((pos = line.find(delimiter)) != std::string::npos) 
        {
            seq = line.substr(0, pos);
            line.erase(0, pos + 1);
            names.push_back(seq);
        }
        seq = line;
        names.push_back(seq);
    }
    nameFileStream.close();
    if(names.empty())
    {
        std::cout << "ERROR: Could not parse any names for guides! Check the guide name file.\n";
        exit(EXIT_FAILURE);
    }

    if(names.size() != seqs.size())
    {
        std::cout << "ERROR: The number of sequences and names for guide reads does not match. Check files for guide sequences and names.\n";
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < names.size(); ++i)
    {
        map.insert(std::make_pair(seqs.at(i), names.at(i)));
    }

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
    double umiThreshold = -1;
    bool scClassConstraint = true;
    bool umiRemoval = true;

    //data for protein(ab) and treatment information
    std::string abFile; 
    int abIdx;
    std::string treatmentFile;
    int treatmentIdx = INT_MAX;
    std::vector<std::string> abBarcodes;
    std::vector<std::string> treatmentBarcodes;

    //data for class information (given e.g. by guide RNA)
    std::string classSeqFile;
    std::string classNameFile;
    std::string guideReadsFile;

    //only used in case we have no CombinatiorialIndexing
    //but a single UMI-like SingleCell ID
    std::string umiSingleCellIdx;

    if(!parse_arguments(argv, argc, inFile, outFile, thread, barcodeFile, barcodeIndices, 
                        umiMismatches, abFile, abIdx, treatmentFile, treatmentIdx,
                        classSeqFile, classNameFile, umiThreshold, scClassConstraint, 
                        guideReadsFile, umiRemoval, umiSingleCellIdx))
    {
        exit(EXIT_FAILURE);
    }

    //make sure we have ETHER a barcode list for CombinatorialIndexing (barcodeIndices)
    // OR a single UMI-like single cell Idx
    if(barcodeIndices == "")
    {
        assert( (umiSingleCellIdx != "") && "We can have ETHER the -c OR the -s flag for CI-barcodes OR a single UMI-like single cell sequence.");
    }
    else
    {
        assert( (umiSingleCellIdx == "") &&  "We can have ETHER the -c OR the -s flag for CI-barcodes OR a single UMI-like single cell sequence.");
    }
    
    //generate the dictionary of barcode alternatives to idx
    NBarcodeInformation barcodeIdData;
    generateBarcodeDicts(barcodeFile, barcodeIndices, barcodeIdData, abBarcodes, abIdx, 
                        umiSingleCellIdx, &treatmentBarcodes, treatmentIdx);
    BarcodeProcessingHandler dataParser(barcodeIdData);
    if(umiThreshold != -1){dataParser.setUmiFilterThreshold(umiThreshold);}
    dataParser.setScClassConstaint(scClassConstraint);
    dataParser.setumiRemoval(umiRemoval);

    //generate dictionaries to map sequences to the real names of Protein/ treatment/ etc...
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
    if(!classSeqFile.empty()) //could take any of both files: generateClassDict checks for both files beeing same length
    {
        std::unordered_map<std::string, std::string > map = generateClassDict(classSeqFile, classNameFile);
        dataParser.addClassData(map);
    }

    //parse the file of demultiplexed barcodes and
    //add all the data to the Unprocessed Demultiplexed Data (stored in rawData)
    // (AB, treatment already are mapped to their real names, scID is a concatenation of numbers for each barcode in
    //each abrcoding round, seperated by a dot)
    if(guideReadsFile != "")
    {
        dataParser.parse_ab_and_guide_file(inFile, guideReadsFile, thread);
    }
    else
    {
        dataParser.parse_combined_file(inFile, thread);
    }
    //further process the data (correct UMIs, collapse same UMIs, etc.)
    dataParser.processBarcodeMapping(umiMismatches, thread);
    dataParser.writeLog(outFile);
    dataParser.writeAbCountsPerSc(outFile);

    return(EXIT_SUCCESS);
}
