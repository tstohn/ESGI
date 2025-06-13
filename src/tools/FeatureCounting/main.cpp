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
                     std::string& barcodeDir, std::string& barcodeIndices, 
                     std::string& umiIdx, int& umiMismatches,
                     std::string& abFile, int& featureIdx, std::string& treatmentFile, int& treatmentIdx,
                     double& umiThreshold, bool& umiRemoval,  bool& scIdString)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input,i", value<std::string>(&inFile)->required(), "input file of demultiplexed reads for ABs in Single cells. (input must be a tsv file, it can be gzipped)")
            ("output,o", value<std::string>(&outFile)->required(), "output file with all split barcodes")

            ("barcodeDir,d", value<std::string>(&(barcodeDir)), " path to a directory which must contain all the barcode files (for variable barcodes). When running <demultiplex> we \
            provided files for variable barcodes in the pattern-file, these files are now in the header of the output of <demultiplex>, but we still need to access those files again to assign features (e.g., protein names) or single-cell IDs to the barcode.")
            
            ("antibodyList,a", value<std::string>(&(abFile)), "file with a list of all feature names (e.g., protein names), should be in same order as the feature-barcodes in the barcode file.")
            ("featureIndex,x", value<int>(&featureIdx)->required(), "Index used for feature counting (e.g., index of the protein barcode). This is the index of the column that should be used for features (0 indexed)")
            
            ("groupList,g", value<std::string>(&(treatmentFile)), "file with a list of all single-cell group assignments (e.g.treatments in specific wells). This is just a file with group names, comma seperated and should be in same order as the specific barcodes for the grouping barcode. \
            If this argument is given, you must also add the index of barcodes used for grouping with <-y>. \
            E.g., imagine we have a barcode file like : ACGT,TACG,CCCG. And the barcodes also define different treatment conditions \
            then we can provide a grouping file -g groupingFile.txt with groupingFile.txt: untreated, treated_time1, treated_time2")
            ("groupingIndex,y", value<int>(&treatmentIdx), "Index used to group cells (e.g. by treatment). This is the barcode used to assign groups that have to be given in <-g>. This is the x-th barcode from the barcodeFile (0 indexed).")

            ("singleCellIndices,c", value<std::string>(&(barcodeIndices))->default_value(""), "comma seperated list of indexes, that are used for \
            single-cell assignment (e.g., for combinatorial indexing 0,5,3. If cells have a single barcode it can be, e.g., only 0). \
            These indices are the indices in the barcode-pattern file (0 indexed). E.g., a pattern like ABPATTERN:[ACGT][5X][scFile1.txt][ACGT][scFile2.txt] could have the single-cell ids -c 2,4.")

            ("umiIndex,u", value<std::string>(&umiIdx)->default_value(""), "list of indices used as unique molecular identifier (UMI). This can be several columns (0 indexed). \
            In the example pattern ABPATTERN:[ACGT][5X][scFile1.txt][ACGT][scFile2.txt] this parameter could be -u 1.\
            If this parameter is not given all columns with an X from the pattern-input-file (e.g., [10X]) are used as UMI. But you can also just provide one pattern to be used as UMI in case several are present.")
            ("mismatches,m", value<int>(&umiMismatches)->default_value(1), "number of allowed mismatches in a UMI (all UMIs are aligned to one another and collapsed if possible). If there are several UMI-barcodes in one sequence\
            the sequences are concatenated and the whole sequence is aligned to other UMI-seuqences by this here provided mismatch number.")
            ("umiThreshold,f", value<double>(&umiThreshold)->default_value(0.0), "threshold for filtering UMIs. E.g. if set to 0.9 we only retain reads of a UMI, if more \
            than 90percent of them have the same SC-AB combination. All other reads are deleted. Default is zero. (You can keep it at 0 if UMIs should not be collapsed).")
            ("umiRemoval,z", value<bool>(&umiRemoval)->default_value(true), "Set to false if UMIs should NOT be collapsed. By default UMIs are collapsed.")
            ("scIdAsString,s", value<bool>(&scIdString)->default_value(false), "Stores the single-cell ID not as an id for the barcode (e.g., 1.45.0), but as the actual string (e.g., ATCG.ACTGT.GCGC).")

            ("thread,t", value<int>(&threats)->default_value(5), "number of threads")
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
std::unordered_map<std::string, std::string > generateProteinDict(std::string abFile, 
                                                                  const std::vector<std::string>& abBarcodes)
{
    std::unordered_map<std::string, std::string > map;
    std::vector<std::string> proteinNames;

    std::ifstream abFileStream(abFile);
    if (!abFileStream.is_open()) 
    {
        std::cerr << "Error: Failed to open antibody file" << abFile << ". Please double check if the file exists.\n";
        exit(EXIT_FAILURE);
    }

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
    for(size_t i = 0; i < abBarcodes.size(); ++i)
    {
        map.insert(std::make_pair(abBarcodes.at(i), proteinNames.at(i)));
    }
    abFileStream.close();

    return map;
}

// generate a dictionary to map sequences to treatments
std::unordered_map<std::string, std::string > generateTreatmentDict(std::string treatmentFile,
                                                                    const std::vector<std::string>& treatmentBarcodes)
{
    std::unordered_map<std::string, std::string > map;
    std::vector<std::string> treatmentNames;

    std::ifstream treatmentFileStream(treatmentFile);
    if (!treatmentFileStream.is_open()) 
    {
        std::cerr << "Error: Failed to open cell-grouping file" << treatmentFile << ". Please double check if the file exists.\n";
        exit(EXIT_FAILURE);
    }

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
    for(size_t i = 0; i < treatmentBarcodes.size(); ++i)
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
    std::string barcodeDir;
    std::string barcodeIndices;
    int thread;
    int umiMismatches;
    double umiThreshold = -1;
    bool umiRemoval = true;
    bool scIdAsString = false;

    //data for protein(ab) and treatment information
    std::string abFile; 
    int featureIdx;
    std::string treatmentFile;
    int treatmentIdx = -1;
    std::vector<std::string> abBarcodes;
    std::vector<std::string> treatmentBarcodes;
    std::string umiIdx;

    if(!parse_arguments(argv, argc, inFile, outFile, thread, 
                        barcodeDir, barcodeIndices, umiIdx, umiMismatches, 
                        abFile, featureIdx, treatmentFile, treatmentIdx,
                        umiThreshold, umiRemoval, scIdAsString))
    {
        exit(EXIT_FAILURE);
    }
    
    //generate the dictionary of barcode alternatives to idx
    BarcodeInformation barcodeIdData;

    //get the first line of headers from input file
    std::ifstream file;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    bool gz = isGzipped(inFile);
    std::istream* instream = openFile(inFile, file, inbuf, gz);
    if (!instream)
    {
        std::cerr << "Error reading input file or file is empty! Please double check if the file exists:" << inFile << std::endl;
        exit(EXIT_FAILURE);
    };

    std::string firstLine;
    if (std::getline(*instream, firstLine)) 
    {  
        bool parseAbBarcodes = true;
        if(abFile.empty()){parseAbBarcodes = false;}
        generateBarcodeDicts(firstLine, barcodeDir, barcodeIndices, barcodeIdData, abBarcodes, parseAbBarcodes, featureIdx, &treatmentBarcodes, treatmentIdx, umiIdx, umiMismatches);
    } 
    else
    {
        std::cerr << "Error reading header line of input file:" << inFile << std::endl;
        exit(EXIT_FAILURE);
    }
    //clean data if necessary
    if (instream != &file) delete instream;
    instream = nullptr;

    BarcodeProcessingHandler dataParser(barcodeIdData);
    if(umiThreshold != -1){dataParser.setUmiFilterThreshold(umiThreshold);}
    dataParser.setumiRemoval(umiRemoval);
    dataParser.setSingleCellIdStyle(scIdAsString);

    //generate dictionaries to map sequences to the real names of Protein/ treatment/ etc...
    std::unordered_map<std::string, std::string > featureMap;
    if(!abFile.empty())
    {
        featureMap = generateProteinDict(abFile, abBarcodes);
    }

    //featureMap is empty if the feature name should stay as they are (no mapping of e.g. barcodes to proteins)
    dataParser.addProteinData(featureMap);
    if(!treatmentFile.empty() && treatmentIdx != -1)
    {
        std::unordered_map<std::string, std::string > treatmentMap = generateTreatmentDict(treatmentFile, treatmentBarcodes);
        dataParser.addTreatmentData(treatmentMap);
    }

    //parse the file of demultiplexed barcodes and
    //add all the data to the Unprocessed Demultiplexed Data (stored in rawData)
    // (AB, treatment already are mapped to their real names, scID is a concatenation of numbers for each barcode in
    //each abrcoding round, seperated by a dot)
    dataParser.parse_barcode_file(inFile);

    //further process the data (correct UMIs, collapse same UMIs, etc.)
    dataParser.processBarcodeMapping(thread);
    dataParser.writeLog(outFile);
    dataParser.writeAbCountsPerSc(outFile);

    return(EXIT_SUCCESS);
}
