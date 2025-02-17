#include "DemultiplexedLinesWriter.hpp"

//initialize the statistics file/ lines that could not be mapped
//FILES: mismatches per barcode / mismatches per barcodePattern/ failedLines
void initialize_additional_output(std::string output)
{
    //remove output
    std::string barcodeMismatches = output + "/BarcodeMismatches.txt";
    std::string patternMismatches = output + "/PatternMismatches.txt";
    std::string failedLines = output + "/failedLines.txt";

    // remove outputfile if it exists
    std::remove(barcodeMismatches.c_str());
    std::remove(patternMismatches.c_str());
    std::remove(failedLines.c_str());
}

//function to replace function below: initialize_output

//this function only initialized the output files that are needed for a specific pattern
//like fastq, demultiplexed-read files
void initialize_output_for_pattern(std::string output, const BarcodePatternPtr pattern)
{
    //txt-file with barcodes
    std::string barcodeOutput = output + "/" + pattern->patternName + ".txt";
    std::remove(barcodeOutput.c_str());

    //fastq-file with the RNA sequence
    std::string rnaOutput;
    if(pattern->containsDNA)
    {
        rnaOutput = output + "/" + pattern->patternName + ".fastq";
        std::remove(rnaOutput.c_str());
    }

    //write header line for barcode file: e.g.: [ACGGCATG][BC1.txt][15X]
    std::ofstream barcodeOutputStream;   
    barcodeOutputStream.open(barcodeOutput, std::ofstream::app);

    for(int bidx = 0; bidx < (pattern->barcodePattern)->size(); ++bidx)
    {

        BarcodePtr bptr = (pattern->barcodePattern)->at(bidx);
        //stop and DNA pattern should not be written
        if( (bptr->name != "*") && (bptr->name != "DNA"))
        {
            barcodeOutputStream << bptr->name;
            if( bidx != ((pattern->barcodePattern)->size()-1))
            {
                barcodeOutputStream << "\t";
            }
        }
    }
    barcodeOutputStream << "\n";
    barcodeOutputStream.close();

}




/// creates new files for failed lines, mapped barcodes (and writes header), statistics
void initialize_output(std::string output, const std::vector<std::pair<std::string, char> > patterns, 
                       std::string& guideNameTage, bool initializeGuideFile = false, bool guideFileHasUmi = false)
{
    //remove output
    std::string outputStats;
    std::string outputMapped;
    std::string outputFailed;
    std::string outputGuide;

    std::size_t found = output.find_last_of("/");
    if(found == std::string::npos)
    {
        outputMapped = "Demultiplexed_" + output;
        outputStats = "StatsMismatches_" + output;
        outputFailed = "FailedLines_" + output;
    }
    else
    {
        outputMapped = output.substr(0,found) + "/" + "Demultiplexed_" + output.substr(found+1);
        outputStats = output.substr(0,found) + "/" + "StatsMismatches_" + output.substr(found+1);
        outputFailed = output.substr(0,found) + "/" + "FailedLines_" + output.substr(found+1);
    }
    // remove outputfile if it exists
    std::remove(outputMapped.c_str());
    std::remove(outputStats.c_str());
    std::remove(outputFailed.c_str());
    std::remove(outputGuide.c_str());

    //write header line for AB file
    std::ofstream outputFile;   
    outputFile.open (outputMapped, std::ofstream::app);
    for(int i =0; i < patterns.size(); ++i)
    {
        //stop pattern should not be written
        if(patterns.at(i).second != 's')
        {
            outputFile << patterns.at(i).first;
            if( i!=(patterns.size() - 1))
            {
            outputFile << "\t";
            }
        }
    }
    outputFile << "\n";
    outputFile.close();

    //write header line for guide file
    if(initializeGuideFile)
    {
        std::ofstream outputFile;   
        outputFile.open (outputGuide, std::ofstream::app);
        for(int i =0; i < patterns.size(); ++i)
        {
            if( (patterns.at(i).second != 'w') || (guideFileHasUmi))
            {
                outputFile << patterns.at(i).first;
                if( i!=(patterns.size() - 1) )
                {
                    outputFile << "\t";
                }
            }
        }
        outputFile << "\n";
        outputFile.close();
    }
}

/// write mismatches per barcode to file
void write_stats(const input& input, const std::map<std::string, std::vector<int> >& statsMismatchDict)
{
    std::string output = input.outPath;
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");
    if(found == std::string::npos)
    {
        output = "StatsMismatches_" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "StatsMismatches_" + output.substr(found+1);
    }
    std::remove(output.c_str());
    outputFile.open (output, std::ofstream::app);

    for(std::pair<std::string, std::vector<int> > mismatchDictEntry : statsMismatchDict)
    {
        outputFile << mismatchDictEntry.first << "\t";
        for(int mismatchCountIdx = 0; mismatchCountIdx < mismatchDictEntry.second.size(); ++mismatchCountIdx)
        {
            outputFile << mismatchDictEntry.second.at(mismatchCountIdx);
            if(mismatchCountIdx!=mismatchDictEntry.second.size()-1){outputFile << "\t";}
        }
        outputFile << "\n";
    }

    outputFile.close();
}

/// write failed lines into a txt file
void write_failed_line(const input& input, std::pair<const std::string&, const std::string&> failedLine)
{
    if(failedLine.second.empty())
    {
        std::string output = input.outPath;
        std::ofstream outputFile;

        //write real sequences that map to barcodes
        std::size_t found = output.find_last_of("/");
        std::string outputFail;
        if(found == std::string::npos)
        {
            outputFail = "FailedLines_" + output;
        }
        else
        {
            outputFail = output.substr(0,found) + "/" + "FailedLines_" + output.substr(found+1);
        }
        outputFile.open (outputFail, std::ofstream::app);
        
        outputFile << failedLine.first << "\n";
        
        outputFile.close();
    }
    else
    {
        //write forward reads
        std::string output = input.outPath;
        std::ofstream outputFile;

        //write real sequences that map to barcodes
        std::size_t found = output.find_last_of("/");
        std::string outputFail;
        if(found == std::string::npos)
        {
            outputFail = "FailedLines_1_" + output;
        }
        else
        {
            outputFail = output.substr(0,found) + "/" + "FailedLines_1_" + output.substr(found+1);
        }
        outputFile.open (outputFail, std::ofstream::app);
        
        outputFile << failedLine.first << "\n";
        
        outputFile.close();

        //write reverse reads
        //write real sequences that map to barcodes
        if(found == std::string::npos)
        {
            outputFail = "FailedLines_2_" + output;
        }
        else
        {
            outputFail = output.substr(0,found) + "/" + "FailedLines_2_" + output.substr(found+1);
        }
        outputFile.open (outputFail, std::ofstream::app);
        
        outputFile << failedLine.second << "\n";
        
        outputFile.close();
    }

}

/// write mapped barcodes to a tab separated file
void write_file(const input& input, BarcodeMappingVector barcodes, std::string nameTag = "")
{
    std::string output = input.outPath;
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");

    //write the barcodes we mapped
    if(found == std::string::npos)
    {
        output = "Demultiplexed_" + nameTag + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "Demultiplexed_" + nameTag + output.substr(found+1);
    }
    outputFile.open (output, std::ofstream::app);
    for(int i = 0; i < barcodes.size(); ++i)
    {
        for(int j = 0; j < barcodes.at(i).size(); ++j)
        {
            outputFile << barcodes.at(i).at(j);
            if(j!=barcodes.at(i).size()-1){outputFile << "\t";}
        }
        outputFile << "\n";
    }
    outputFile.close();
}

/// calls output initializer functions and gets the barcode mapping structure from Mapping object, since this will the header of the output file
// create backbone files for barcoding patterns that will be mapped: e.g.: FASTQ for RNA, txt with heads for barcode-files for CI, spatial, other stuff
//the file will be anmed after pattern name
template <typename MappingPolicy, typename FilePolicy>
void DemultiplexedLinesWriter<MappingPolicy, FilePolicy>::initialize_output_files(const input& input)
{
    //TO DO
    //parse through the barcodePatterns, make file of pattern name
    for(const BarcodePatternPtr barcodePattern : *this->barcodePatternList)
    {
        initialize_output_for_pattern(input.outPath, barcodePattern);
    }

    //create universal output files
    // 2 statistics files: mismatches per pattern, and mismatches in barcodes
    // file with failed lines
    initialize_additional_output(input.outPath);
}


/**
* @brief function wrapping the demultiplex_read function of the Mapping class to increment and decrement a counter of elements
* in the current queue. Used to allow processing of only a few lines a time instead of writing all into RAM.
**/
template <typename MappingPolicy, typename FilePolicy>
void DemultiplexedLinesWriter<MappingPolicy, FilePolicy>::demultiplex_wrapper(std::pair<const std::string&, const std::string&> line,
                                                            const input& input,
                                                            std::atomic<unsigned long long>& lineCount,
                                                            const unsigned long long& totalReadCount,
                                                            std::atomic<long long int>& elementsInQueue)
{

    //TO DO:
    //here iterate through the list of possible barcodePatterns and map one after the other
    //so far was done in barcdoeMapping - so remove from there and hadnle the barcodePattern options here

    //firstly try mapping an AB read
    bool result = this->demultiplex_read(line, input, lineCount, totalReadCount, false);

    /*if(!result && input.guideFile != "")
    {
        //run again this time mapping guide reads
        result = this->demultiplex_read(line, input, lineCount, totalReadCount, true);
    }
*/

    if(!result && input.writeFailedLines)
    {
        //write failed line to file
        write_failed_line(input, line);
    }
    --elementsInQueue;


}

/// overwritten run_mapping function to allow processing of only a subset of fastq lines at a time
template <typename MappingPolicy, typename FilePolicy>
void DemultiplexedLinesWriter<MappingPolicy, FilePolicy>::run_mapping(const input& input)
{
    std::cout << "START DEMULTIPLEXING\n";

    //generate a pool of threads
    boost::asio::thread_pool pool(input.threads); //create thread pool

    //read line by line and add to thread pool
    this->FilePolicy::init_file(input.inFile, input.reverseFile);
    std::pair<std::string, std::string> line;
    std::atomic<unsigned long long> lineCount = 0; //using atomic<int> as thread safe read count
    std::atomic<long long int> elementsInQueue = 0;
    unsigned long long totalReadCount = FilePolicy::get_read_number();

    // have every thread keep its own copy of *this*object , to write barcodes without locking
    //then in the end have every thread write its mappings to a file
    while(FilePolicy::get_next_line(line))
    {
        //wait to enqueue new elements in case we have a maximum bucket size
        if(input.fastqReadBucketSize>0)
        {
            while(input.fastqReadBucketSize <= elementsInQueue){}
        }
        //increase job count and push the job in the queue
        ++elementsInQueue;
        boost::asio::post(pool, std::bind(&DemultiplexedLinesWriter::demultiplex_wrapper, this, line, input, std::ref(lineCount), totalReadCount, std::ref(elementsInQueue)));
    }

    //iterate through the barcode map and let threads 
    pool.join();
    printProgress(1); std::cout << "\n"; // end the progress bar
    if(totalReadCount != ULLONG_MAX)
    {
        std::cout << "=>\tPERFECT MATCHES: " << std::to_string((unsigned long long)(100*(this->get_perfect_matches())/(double)totalReadCount)) 
                << "% | MODERATE MATCHES: " << std::to_string((unsigned long long)(100*(this->get_moderat_matches())/(double)totalReadCount))
                << "% | MISMATCHES: " << std::to_string((unsigned long long)(100*(this->get_failed_matches())/(double)totalReadCount)) << "%\n";
    }

    FilePolicy::close_file();
}

/**
* @brief overwritten run_mapping function of Mapping class to allow processing of only a subset of fastq lines at a time
* and to store all output results that we want to safe (e.g. failed lines, statistics)
**/
template <typename MappingPolicy, typename FilePolicy>
void DemultiplexedLinesWriter<MappingPolicy, FilePolicy>::run(const input& input)
{
    //from the basic information within patterns generate a more complex barcodePattern object
    //which stores for each pattern all possible barcodes, number of mismatches etc.
    this->generate_barcode_patterns(input);

    //create output files and write headers for demultiplexed barcodes
    initialize_output_files(input);

    //create empty dict for mismatches per barcode
    if(input.writeStats)
    {
        this->initializeStats();
    }

    //run mapping
    this->run_mapping(input);

    //write the barcodes, failed lines, statistics (mismatches per barcode)
    //TODO:
    //iterate over the different result maps for the various barcodePatterns

    write_file(input, this->get_demultiplexed_ab_reads());
    write_stats(input, this->get_mismatch_dict());
}

template class DemultiplexedLinesWriter<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy>;
template class DemultiplexedLinesWriter<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy>;
template class DemultiplexedLinesWriter<MapEachBarcodeSequentiallyPolicyPairwise, ExtractLinesFromFastqFilePolicyPairedEnd>;