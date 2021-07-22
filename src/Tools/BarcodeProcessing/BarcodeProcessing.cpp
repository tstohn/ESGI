#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

#include "AbData.hpp"

using namespace boost::program_options;

/**
 * HOW TO:
 * PREREQUISITS:
 *  THE BARCODE FILE TO ANALYSE MUST HAVE A HEADER WITH THE BARCODEPATTERNS AS IN FASTQPARSER TOOL, THERE MUST BE A UMI SEQU< AND A AB SEQ, which MUST BE THE ONLY
 *  VARYING BARCODE SEQUENCE EXCEPT THE CI BARCODES
 * 
 * - generates two files: 
 *          - firstly collapse all UMI_barcodes and have a file: UMI_idx, Ab_idx, cell_idx
 *          - secondly collapse all UMI_idxs: Ab_idx, Ab_count, cell_idx
 * 
 *  - get all unique UMIs
 *  - cluster UMIs that r similar to the same group of UMIs
 * => make a file where UMI column is the 'expected' UMI, from this file we can calculate the saturation...
 * 
 *  - delete doublicate reads (based on UMI and BC 1-4)
 *  - for the same cell (based on BC 1-4) count the number of UMIs per AB_barcode and write a new column AB_count,
 *    so that we end up with a file Ab_id, sample_id, Ab_count
 * => generate the file format for our ScRNA_seq_Normalization pipeline
 * 
 * DISTINGUISH THOSE TWO BY UMI-DATA and BARCODE-DATA
 * 
 * PARAMETER FOR ABRCODE POSITIONS:
 *  give a list of indices for the CIBarcodes within the file, which hold all alternatives for the varying indices (the real position indide the fastqRead is then parsed
 *  by checking for NNNN sequences), the other one left NNN sequence is the AB sequence and the XXX is the UMI sequence
 * */

//statistics of the UMI occurences
struct StatsUmi
{
    
};

bool parse_arguments(char** argv, int argc, std::string& inFile,  std::string& outFile, int& threats, std::string& barcodeFile, std::string& barcodeIndices)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input,i", value<std::string>(&inFile)->required(), "directory of files or single file in fastq(.gz) format")
            ("output,o", value<std::string>(&outFile)->required(), "output file with all split barcodes")
            ("barcodeList,b", value<std::string>(&(barcodeFile)), "file with a list of all allowed well barcodes (comma seperated barcodes across several rows)\
            the row refers to the correponding bracket enclosed sequence substring. E.g. for two bracket enclosed substrings in out sequence a possible list could be:\
            AGCTTCGAG,ACGTTCAGG\nACGTCTAGACT,ATCGGCATACG,ATCGCGATC,ATCGCGCATAC. This can be the same list as it was for FastqParser.")
            
            ("CombinatorialIndexingBarcodeIndices,c", value<std::string>(&(barcodeIndices))->required(), "comma seperated list of indexes, that are used during \
            combinatorial indexing and should distinguish a unique cell. Be aware that this is the index of the line inside the barcodeList file (see above). \
            This file ONLY includes lines for the varying sequences (except UMI). Therefore the index is not the same as the position in the whole sequence \
            if constant or UMI-seq are present. Index starts with zero.")

            ("thread,t", value<int>(&threats)->default_value(5), "number of threads")
            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";
            std::cout << "EXAMPLE CALL:\n ./bin/processing -i <inFile>\n";
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

void generateBarcodeDicts(std::string barcodeFile, std::string barcodeIndices, CIBarcode& barcodeIdData)
{
    //parse barcode file
    std::vector<std::vector<std::string> > barcodeList;
    std::ifstream barcodeFileStream(barcodeFile);
    for(std::string line; std::getline(barcodeFileStream, line);)
    {
        std::string delimiter = ",";
        std::string seq;
        size_t pos = 0;
        std::vector<std::string> seqVector;
        while ((pos = line.find(delimiter)) != std::string::npos) 
        {
            seq = line.substr(0, pos);
            line.erase(0, pos + 1);
            for (char const &c: seq) {
                if(!(c=='A' | c=='T' | c=='G' |c=='C' |
                        c=='a' | c=='t' | c=='g' | c=='c'))
                        {
                        std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C)\n";
                        if(c==' ' | c=='\t' | c=='\n')
                        {
                            std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
                        }
                        exit(1);
                        }
            }
            seqVector.push_back(seq);
        }
        seq = line;
        for (char const &c: seq) {
            if(!(c=='A' | c=='T' | c=='G' |c=='C' |
                    c=='a' | c=='t' | c=='g' | c=='c'))
                    {
                    std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C)\n";
                    if(c==' ' | c=='\t' | c=='\n')
                    {
                        std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
                    }
                    exit(1);
                    }
        }
        seqVector.push_back(seq);
        barcodeList.push_back(seqVector);
        seqVector.clear();
    }
    barcodeFileStream.close();

    //parse the indices of CIbarcodes from line
    std::stringstream ss;
    ss.str(barcodeIndices);
    while(ss.good())
    {
        std::string substr;
        getline(ss, substr, ',' );
        barcodeIdData.ciBarcodeIndices.push_back(stoi(substr));
    }

    //vector of barodeList stores all possible barcodes at each index
    //generate a struct with a dictionary maping a barcode to an index
    //for every CI Barcode in sequence
    for(const int& i : barcodeIdData.ciBarcodeIndices)
    {
        int barcodeCount = 0;
        std::unordered_map<std::string, int> barcodeMap;
        //for all options of this barcode
        for(const std::string& barcodeEntry : barcodeList.at(i))
        {
            barcodeMap.insert(std::pair<std::string, int>(barcodeEntry,barcodeCount));
            ++barcodeCount;
        }
        barcodeIdData.barcodeIdDict.push_back(barcodeMap);
    }

}


int main(int argc, char** argv)
{
    std::string inFile;
    std::string outFile;
    std::string barcodeFile;
    std::string barcodeIndices;
    int thread;
    parse_arguments(argv, argc, inFile, outFile, thread, barcodeFile, barcodeIndices);
    
    //generate the dictionary of barcode alternatives to idx
    CIBarcode barcodeIdData;
    generateBarcodeDicts(barcodeFile, barcodeIndices, barcodeIdData);
    
    UmiDataParser dataParser(barcodeIdData);
    dataParser.parseFile(inFile, thread);

    //dataParser.writeStats(outFile);
    dataParser.writeUmiCorrectedData(outFile);
    //AbData abData(data);
    //abData.writeFile(outFile);

    return(EXIT_SUCCESS);
}