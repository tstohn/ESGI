#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <zlib.h>
#include <regex>
#include <thread>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

#include "UmiDataParser.hpp"
#include "mapping.hpp"

using namespace boost::program_options;

/*
    1.) UMI check
         if same UMI and NOT same AB-SC: check if same AB, check if same BC 4 1 2 3

    2.) Barcode Pattern Check



*/

struct umiRead
{
    const char* sc;
    const char* ab;
};

bool parse_arguments(char** argv, int argc, std::string& inputFails, std::string& inputBarcodeMapping, std::string& output, 
                     int& abIdx, int& treatmentIdx, int& threads, std::string& barcodeFile, std::string& barcodeIndices,
                     std::string& seqpattern, std::string& mismatches,std::string& correctredUmiFile)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("output,o", value<std::string>(&(output))->required(), "output file with all split barcodes")
            ("barcodeList,b", value<std::string>(&(barcodeFile)), "file with a list of all allowed well barcodes (comma seperated barcodes across several rows)\
            the row refers to the correponding bracket enclosed sequence substring. E.g. for two bracket enclosed substrings in out sequence a possible list could be:\
            AGCTTCGAG,ACGTTCAGG\nACGTCTAGACT,ATCGGCATACG,ATCGCGATC,ATCGCGCATAC. This can be the same list as it was for FastqParser.")

            //parameters for failed lines
            ("inputFails,f", value<std::string>(&inputFails), "lines of nucleotide sequences that failed the analysis")
            ("sequencePattern,p", value<std::string>(&seqpattern), "pattern for the sequence to match, \
            every substring that should be matched is enclosed with square brackets. N is a wild card match, a substring \
            should not be a combination of wild card chars and constant chars : [AGCTATCACGTAGC][NNNNNN][AGAGCATGCCTTCAG][NNNNNN]")
            ("mismatches,m", value<std::string>(&(mismatches)), "list of mismatches allowed for each bracket enclosed sequence substring. \
            This should be a comma seperated list of numbers for each substring of the sequence enclosed in squared brackets. E.g.: 2,1,2,1,2")

            //parameters for UMI CHECK
            ("inputBarcodes,i", value<std::string>(&inputBarcodeMapping), "lines of tab seperated barcode mappings")
            ("antibodyIndex,x", value<int>(&abIdx), "Index used for antibody distinction.")
            ("CombinatorialIndexingBarcodeIndices,c", value<std::string>(&(barcodeIndices)), "comma seperated list of indexes, that are used during \
            combinatorial indexing and should distinguish a unique cell. Be aware that this is the index of the line inside the barcodeList file (see above). \
            This file ONLY includes lines for the varying sequences (except UMI). Therefore the index is not the same as the position in the whole sequence \
            if constant or UMI-seq are present. Index starts with zero.")
            ("correctredUmiFile,u", value<std::string>(&correctredUmiFile), "output file of barcode mapping with the corrected UMI reads (+ AB and SC ids)")

            ("thread,t", value<int>(&threads)->default_value(1), "number of threads")

            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";

            std::cout << "###########################################\n";
            std::cout << "EXAMPLE CALL:\n ./bin/parser -i ./inFile -o ./outFile -p [AGTCAGTC][NNNN] -b ./barcodeFile.txt -m 2,1 -t 5\n";
            std::cout << "Calling the Tool, mapping each line to a pattern, that starts with a constant sequence of <AGTCAGTC> in which up to two mismatches\n\
            are allowed. After that the reads should contain a four base long sequence with a maximum of one mismatch. All sequences that are allowed are in the text file barcodeFile.txt.\n\
            This file should have in its first row, all comma seperated sequences that could map, and they should be all only four bases long, allowed chars are only A,C,G,T.\n";
            std::cout << "###########################################\n";

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

void analyse_failed_lines(input& input)
{

    BarcodeMappingVector barcodes;
    BarcodeMappingVector realBarcodes;
    std::vector<std::string> fastqLines;
    std::vector<std::pair<std::string, char> > patterns; // vector of all string patterns, 
                                                        //second entry is c=constant, v=varying, w=wildcard
    BarcodePatternVectorPtr barcodePatterns = generate_barcode_patterns(input, patterns);

    //initialize stats
    fastqStats fastqStatsFinal;
    initializeStats(fastqStatsFinal, barcodePatterns);
    initializeOutput(input.outFile, patterns);

    //read lines into memory
    std::ifstream fileStream(input.inFile);
    int totalReads = std::count(std::istreambuf_iterator<char>(fileStream), std::istreambuf_iterator<char>(), '\n');
    fileStream.clear();
    fileStream.seekg(0);
    bool readsLeft = true;
    int totalCurrentReadNum = 0;
    fastqStats emptyStats = fastqStatsFinal;
    while(readsLeft)
    {
        //push a batch of reads into a temporary vector
        fastqLines.clear();
        int fastqReadThreashold = input.fastqReadBucketSize;
        int numFastqReads = 0;
        while(numFastqReads<fastqReadThreashold)
        {
            //get next fastq line
            std::string line;
            if(! std::getline(fileStream, line))
            {
                readsLeft = false;
                numFastqReads = fastqReadThreashold;
                continue;
            }
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());

            fastqLines.push_back(line);
            ++numFastqReads;
            ++totalCurrentReadNum;

            if((totalReads >= 100) && (totalCurrentReadNum % (totalReads / 100) == 0))
            {
                double perc = totalCurrentReadNum/ (double)totalReads;
                //printProgress(perc);
            }

        }

        split_barcodes_in_subset(input, barcodes, realBarcodes, barcodePatterns, fastqStatsFinal, patterns, fastqLines, emptyStats);

    }

}

void analyse_same_umis(std::string& inputBarcodeMapping, std::string& output, int& abIdx, const int& threads,
                    const std::string& barcodeFile, const std::string& barcodeIndices)
{
    //data for protein(ab) and treatment information
    std::string abFile; 
    std::string treatmentFile;
    std::vector<std::string> abBarcodes;

    //generate the dictionary of barcode alternatives to idx
    CIBarcode barcodeIdData;

    generateBarcodeDicts(barcodeFile, barcodeIndices, barcodeIdData, abBarcodes, abIdx);

    UmiDataParser dataParser(barcodeIdData);
    //parse the information

    dataParser.parseFile(inputBarcodeMapping, threads);

    //analyse all UMIS
    dataParser.extended_umi_quality_check(threads, output);
}

void analyse_corrected_umis(const std::string& correctredUmiFile, const std::string& output)
{
    //remove previous file and open for reading as append later
    std::ifstream barcodeFileStream(correctredUmiFile);
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");
    std::string umiOutput = output;
    if(found == std::string::npos)
    {
        umiOutput = "QC_CorrectedUmiDist_" + output;
    }
    else
    {
        umiOutput = output.substr(0,found) + "/" + "QC_CorrectedUmiDist_" + output.substr(found+1);
    }
    std::remove(umiOutput.c_str());
    outputFile.open (umiOutput, std::ofstream::app);
    outputFile << "UMI\tREAD_COUNT\tPERCENTAGES\n";
    outputFile.close();

    //dict storing for each Umi all reads
    std::unordered_map<const char*, std::vector<umiRead>, CharHash, CharPtrComparator> umiToDatalines;
    UniqueCharSet uniqueChars;
    int currentLineCount = 0;
    int totalDataSize = totalNumberOfLines(correctredUmiFile);

    std::cout << "Reading Data into Memory\n";
    for(std::string line; std::getline(barcodeFileStream, line);)
    {
        std::string delimiter = "\t";
        std::vector<std::string> splitLine = splitByDelimiter(line, delimiter);
        if(strcmp(splitLine.at(0).c_str(),"UMI") == 0){continue;}
        //order in vector: UMI AB_BARCODE SC_BCID_COMBINATION
        umiRead newRead;
        newRead.ab = uniqueChars.getUniqueChar(splitLine.at(1).c_str());
        newRead.sc = uniqueChars.getUniqueChar(splitLine.at(2).c_str());

        if(umiToDatalines.find(uniqueChars.getUniqueChar(splitLine.at(0).c_str())) == umiToDatalines.end())
        {
            std::vector<umiRead> vec;
            vec.push_back(newRead);
            umiToDatalines.insert(std::make_pair(uniqueChars.getUniqueChar(splitLine.at(0).c_str()), vec));
        }
        //if not add this new umi with this actual position to map
        else
        {
            umiToDatalines[uniqueChars.getUniqueChar(splitLine.at(0).c_str())].push_back(newRead);
        }

        ++currentLineCount;
        if( (totalDataSize >= 100) && (currentLineCount % (totalDataSize / 100) == 0) )
        {
            double perc = currentLineCount/ (double) totalDataSize;
            printProgress(perc);
        }
    }

    //for each umi calculate a vector of percentages
    std::cout << "\nCalculate read distributions for UMI\n";
    currentLineCount = 0;
    totalDataSize = umiToDatalines.size();
    for (std::unordered_map<const char*, std::vector<umiRead>, CharHash, CharPtrComparator>::iterator umiIt = umiToDatalines.begin(); 
         umiIt != umiToDatalines.end(); ++umiIt)
    {
        const char* umiSeq = umiIt->first;
        std::vector<umiRead> reads = umiIt->second;
        unsigned long totalCount = reads.size();
        std::vector<double> percReads;
        while(!reads.empty())
        {
            std::vector<umiRead> newReads;
            umiRead read = reads.at(0);
            unsigned long readCount = 1;
            for(int i =1; i < reads.size(); ++i)
            {

                if( read.ab == reads.at(i).ab && read.sc == reads.at(i).sc)
                {
                    ++readCount;
                }
                else
                {
                    newReads.push_back(reads.at(i));
                }
            }
            reads = newReads;
            percReads.push_back(double(readCount)/totalCount);
        }

        // write UMI, numer of reads and percentages to file
        outputFile.open (umiOutput, std::ofstream::app);
        for(int i = 0; i < percReads.size(); ++i)
        {
                outputFile << umiSeq << "\t" << totalCount << "\t" << percReads.at(i) <<  "\n"; 
        }
        outputFile.close();

        ++currentLineCount;
        if( (totalDataSize >= 100) && (currentLineCount % (totalDataSize / 100) == 0) )
        {
            double perc = currentLineCount/ (double) totalDataSize;
            printProgress(perc);
        }
    }
}

int main(int argc, char** argv)
{
    std::string inputFails;
    std::string inputBarcodeMapping;
    std::string output;

    std::string barcodeFile;
    std::string barcodeIndices;
    std::string seqpattern;
    std::string correctredUmiFile;

    int abIdx;
    int treatmentIdx;
    int threads;
    std::string mismatches;
    if(parse_arguments(argv, argc, inputFails, inputBarcodeMapping, output, abIdx, treatmentIdx, threads, 
                       barcodeFile, barcodeIndices, seqpattern, mismatches, correctredUmiFile))
    {
        //analyse the read number per UMI, distribution of BC combinations, etc
    //  analyse_same_umis(inputBarcodeMapping, output, abIdx, threads, barcodeFile, barcodeIndices);
        //after BarcodeProcessing (UMI MisMatch correction) analyze the occurence of UMIs: number of reads, BC combination percentages
    //  analyse_corrected_umis(correctredUmiFile, output);

        //analze lines that could not be mapped
        input input;
        input.barcodeFile = barcodeFile;
        input.inFile = inputFails;
        input.outFile = output;
        input.patternLine = seqpattern;
        input.mismatchLine = mismatches;
        input.analyseUnmappedPatterns = true;
        input.threads = threads;

        analyse_failed_lines(input);
    }
 
    return EXIT_SUCCESS;
}