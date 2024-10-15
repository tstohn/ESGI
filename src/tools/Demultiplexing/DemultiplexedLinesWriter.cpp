#include "DemultiplexedLinesWriter.hpp"

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
        outputGuide = "Demultiplexed_" + guideNameTage + output;
    }
    else
    {
        outputMapped = output.substr(0,found) + "/" + "Demultiplexed_" + output.substr(found+1);
        outputStats = output.substr(0,found) + "/" + "StatsMismatches_" + output.substr(found+1);
        outputFailed = output.substr(0,found) + "/" + "FailedLines_" + output.substr(found+1);
        outputGuide = output.substr(0,found) + "/" + "Demultiplexed_" + guideNameTage + output.substr(found+1);
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
        outputFile << patterns.at(i).first;
        if( i!=(patterns.size() - 1) )
        {
           outputFile << "\t";
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
    std::string output = input.outFile;
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
        std::string output = input.outFile;
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
        std::string output = input.outFile;
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
    std::string output = input.outFile;
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
template <typename MappingPolicy, typename FilePolicy>
void DemultiplexedLinesWriter<MappingPolicy, FilePolicy>::initialize_output_files(const input& input, 
                                                                                  const std::vector<std::pair<std::string, char> >& patterns,
                                                                                  std::string& guideNameTage)
{
    BarcodePatternVectorPtr barcodePatterns = this->get_barcode_pattern_vector();
    //initialize all output files: write header, delete old files etc.
    bool initializeGuideFile = false, guideFileHasUmi = false;
    if(input.guideFile != ""){initializeGuideFile = true;}
    if(input.guideUMI){guideFileHasUmi = true;}

    initialize_output(input.outFile, patterns, guideNameTage, initializeGuideFile, guideFileHasUmi);
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
    //firstly try mapping an AB read
    bool result = this->demultiplex_read(line, input, lineCount, totalReadCount, false);
    if(!result && input.guideFile != "")
    {
        //run again this time mapping guide reads
        result = this->demultiplex_read(line, input, lineCount, totalReadCount, true);
    }
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
    std::vector<std::pair<std::string, char> > pattern = this->generate_barcode_patterns(input);

    //create output files and write headers for demultiplexed barcodes
    std::string guideNameTage = "guideReads";
    initialize_output_files(input, pattern, guideNameTage);

    //create empty dict for mismatches per barcode
    if(input.writeStats)
    {
        this->initializeStats();
    }

    //run mapping
    this->run_mapping(input);

    //write the barcodes, failed lines, statistics (mismatches per barcode)
    write_file(input, this->get_demultiplexed_ab_reads());
    write_file(input, this->get_demultiplexed_guide_reads(), guideNameTage);
    write_stats(input, this->get_mismatch_dict());
}

template class DemultiplexedLinesWriter<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy>;
template class DemultiplexedLinesWriter<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy>;
template class DemultiplexedLinesWriter<MapEachBarcodeSequentiallyPolicyPairwise, ExtractLinesFromFastqFilePolicyPairedEnd>;