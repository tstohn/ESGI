#include "DemultiplexedLinesWriter.hpp"

//initialize the dictionary of mismatches in each specific barcode
// the vector in the dict has length "mismatches + 2" for each barcode
// one entry for zero mismatches, eveery number from, 1 to mismatches and 
//one for more mismatches than the max allowed number of mismathces
void initializeStats(fastqStats& stats, const BarcodePatternVectorPtr barcodePatterns)
{
    for(BarcodePatternVector::iterator patternItr = barcodePatterns->begin(); 
        patternItr < barcodePatterns->end(); 
        ++patternItr)
    {
        int mismatches = (*patternItr)->mismatches;
        const std::vector<std::string> patterns = (*patternItr)->get_patterns();
        for(const std::string& pattern : patterns)
        {
            std::vector<int> mismatchVector(mismatches + 2, 0);
            stats.mapping_dict.insert(std::make_pair(pattern, mismatchVector));
        }
    }
}

void initializeOutput(std::string output, const std::vector<std::pair<std::string, char> > patterns)
{
    //remove output
    std::string outputStats;
    std::string outputReal;
    std::string outputMapped;
    std::string outputFailed;
    std::size_t found = output.find_last_of("/");
    if(found == std::string::npos)
    {
        outputMapped = "BarcodeMapping_" + output;
        outputStats = "StatsBarcodeMappingErrors_" + output;
        outputReal = "RealBarcodeSequence_" + output;
        outputFailed = "FailedLines_" + output;
    }
    else
    {
        outputMapped = output.substr(0,found) + "/" + "BarcodeMapping_" + output.substr(found+1);
        outputStats = output.substr(0,found) + "/" + "StatsBarcodeMappingErrors_" + output.substr(found+1);
        outputReal = output.substr(0,found) + "/" + "RealBarcodeSequence_" + output.substr(found+1);
        outputFailed = output.substr(0,found) + "/" + "FailedLines_" + output.substr(found+1);
    }
    // remove outputfile if it exists
    std::remove(outputMapped.c_str());
    std::remove(outputStats.c_str());
    std::remove(outputReal.c_str());
    std::remove(outputFailed.c_str());

    //write header line
    std::ofstream outputFile;   
    outputFile.open (outputMapped, std::ofstream::app);
    for(int i =0; i < patterns.size(); ++i)
    {
        outputFile << patterns.at(i).first;
        if( i!=(patterns.size() - 1) )
        {
           outputFile << "\t";
        }
    }
    outputFile << "\n";
    outputFile.close();
}
/*
void write_file(const input& input, BarcodeMappingVector barcodes, BarcodeMappingVector realBarcodes)
{
    std::string output = input.outFile;
    std::ofstream outputFile;

    //write real sequences that map to barcodes
    std::size_t found = output.find_last_of("/");
    if(input.storeRealSequences)
    {
        std::string outputReal;
        if(found == std::string::npos)
        {
            outputReal = "RealBarcodeSequence_" + output;
        }
        else
        {
            outputReal = output.substr(0,found) + "/" + "RealBarcodeSequence_" + output.substr(found+1);
        }
        outputFile.open (outputReal, std::ofstream::app);
        for(int i = 0; i < barcodes.size(); ++i)
        {
            for(int j = 0; j < barcodes.at(i).size(); ++j)
            {
                outputFile << *(barcodes.at(i).at(j));
                if(j!=barcodes.at(i).size()-1){outputFile << "\t";}
            }
            outputFile << "\n";
        }
        outputFile.close();
    }

    //write the barcodes we mapped
    if(found == std::string::npos)
    {
        output = "BarcodeMapping_" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "BarcodeMapping_" + output.substr(found+1);
    }
    outputFile.open (output, std::ofstream::app);
    for(int i = 0; i < realBarcodes.size(); ++i)
    {
        for(int j = 0; j < realBarcodes.at(i).size(); ++j)
        {
            outputFile << *(realBarcodes.at(i).at(j));
            if(j!=barcodes.at(i).size()-1){outputFile << "\t";}
        }
        outputFile << "\n";
    }
    outputFile.close();
}*/

void tmp_quickNdirty_write(const input& input, BarcodeMappingVector barcodes)
{
    std::string output = input.outFile;
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");

    //write the barcodes we mapped
    if(found == std::string::npos)
    {
        output = "BarcodeMapping_" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "BarcodeMapping_" + output.substr(found+1);
    }
    outputFile.open (output, std::ofstream::app);

    for(int i = 0; i < barcodes.size(); ++i)
    {
        for(int j = 0; j < barcodes.at(i).size(); ++j)
        {
            outputFile << (barcodes.at(i).at(j));
            if(j!=barcodes.at(i).size()-1){outputFile << "\t";}
        }
        outputFile << "\n";
    }
    outputFile.close();
}


template <typename MappingPolicy, typename FilePolicy>
void DemultiplexedLinesWriter<MappingPolicy, FilePolicy>::initialize_mapping(const input& input, const std::vector<std::pair<std::string, char> >& patterns)
{
    BarcodePatternVectorPtr barcodePatterns = this->get_barcode_pattern_vector();

    //initialize all output files: write header, delete old files etc.
    initializeOutput(input.outFile, patterns);

    //if stats are written also initialize this file
    if(input.withStats)
    {
        //fastqStats fastqStatsTmp;
        //initializeStats(fastqStatsTmp, barcodePatterns);
        //fastqStatsPtr = std::make_shared<fastqStats>(fastqStatsTmp);
    }

}
/*
template <typename MappingPolicy, typename FilePolicy>
void DemultiplexedLinesWriter<MappingPolicy, FilePolicy>::run_mapping(const input& input)
{
    
    FilePolicy::init_file(input.inFile);

    std::vector<std::string> fastqLines;
    fastqStats emptyStats = *fastqStatsPtr; //an empty statistic with already the right keys for each barcode in the dictionary, used to initialize the stats for each thread

    bool readsLeft = true;
    int totalCurrentReadNum = 0;
    while(readsLeft)
    {
        //push a batch of reads into a temporary vector
        fastqLines.clear();
        int fastqReadThreashold = input.fastqReadBucketSize;
        int numFastqReads = 0;
        while(numFastqReads<fastqReadThreashold && readsLeft)
        {

            std::string line;
            readsLeft = FilePolicy::get_next_line(line);
            if(!readsLeft){continue;}

            fastqLines.push_back(line);
            ++numFastqReads;
            ++totalCurrentReadNum;

            if((FilePolicy::totalReads >= 100) && (totalCurrentReadNum % (FilePolicy::totalReads / 100) == 0))
            {
                double perc = totalCurrentReadNum/ (double)FilePolicy::totalReads;
                printProgress(perc);
            }

        }
        ditribute_jobs_to_threads(input, fastqLines);

    }
    printProgress(1);
    std::cout << "\nMATCHED: " << fastqStatsPtr->perfectMatches << " | MODERATE MATCH: " << fastqStatsPtr->moderateMatches
              << " | MISMATCHED: " << fastqStatsPtr->noMatches << " | Multiplebarcode: " << fastqStatsPtr->multiBarcodeMatch << "\n";

    FilePolicy::close_file();

    writeStats(input.outFile, *fastqStatsPtr);
}*/

template <typename MappingPolicy, typename FilePolicy>
void DemultiplexedLinesWriter<MappingPolicy, FilePolicy>::run(const input& input)
{
    //from the basic information within patterns generate a more complex barcodePattern object
    //which stores for each pattern all possible barcodes, number of mismatches etc.
    std::vector<std::pair<std::string, char> > pattern = this->generate_barcode_patterns(input);

    initialize_mapping(input, pattern);

    this->run_mapping(input);
    tmp_quickNdirty_write(input, this->get_demultiplexed_reads());
}


template class DemultiplexedLinesWriter<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy>;
template class DemultiplexedLinesWriter<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy>;