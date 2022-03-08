#include "MappingAroundLinker.hpp"

/// creates new files for failed lines, mapped barcodes (and writes header), statistics
void initialize_output(std::string output, const std::vector<std::pair<std::string, char> > patterns)
{
    //remove output
    std::string outputMapped;
    std::size_t found = output.find_last_of("/");
    if(found == std::string::npos)
    {
        outputMapped = "DemultiplexedAroundLinker_" + output;
    }
    else
    {
        outputMapped = output.substr(0,found) + "/" + "DemultiplexedAroundLinker_" + output.substr(found+1);
    }
    // remove outputfile if it exists
    std::remove(outputMapped.c_str());

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

/// write mapped barcodes to a tab separated file
void write_file(const input& input, BarcodeMappingVector barcodes)
{
    std::string output = input.outFile;
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");

    //write the barcodes we mapped
    if(found == std::string::npos)
    {
        output = "DemultiplexedAroundLinker_" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "DemultiplexedAroundLinker_" + output.substr(found+1);
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
void MappingAroundLinker<MappingPolicy, FilePolicy>::initialize_output_files(const input& input, const std::vector<std::pair<std::string, char> >& patterns)
{
    BarcodePatternVectorPtr barcodePatterns = this->get_barcode_pattern_vector();
    //initialize all output files: write header, delete old files etc.
    initialize_output(input.outFile, patterns);
}


/**
* @brief function wrapping the demultiplex_read function of the Mapping class to increment and decrement a counter of elements
* in the current queue. Used to allow processing of only a few lines a time instead of writing all into RAM.
**/
template <typename MappingPolicy, typename FilePolicy>
void MappingAroundLinker<MappingPolicy, FilePolicy>::demultiplex_wrapper(std::pair<const std::string&, const std::string&> line,
                                                            const input& input,
                                                            std::atomic<unsigned long long>& lineCount,
                                                            const unsigned long long& totalReadCount,
                                                            std::atomic<long long int>& elementsInQueue)
{
    this->demultiplex_read(line, input, lineCount, totalReadCount);
    --elementsInQueue;
}

/// overwritten run_mapping function to allow processing of only a subset of fastq lines at a time
template <typename MappingPolicy, typename FilePolicy>
void MappingAroundLinker<MappingPolicy, FilePolicy>::run_mapping(const input& input)
{
    std::cout << "START DEMULTIPLEXING OF IMPERFECT BARCODE SEQUENCES\n";

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
        boost::asio::post(pool, std::bind(&MappingAroundLinker::demultiplex_wrapper, this, line, input, std::ref(lineCount), totalReadCount, std::ref(elementsInQueue)));
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
void MappingAroundLinker<MappingPolicy, FilePolicy>::run(const input& input)
{
    //from the basic information within patterns generate a more complex barcodePattern object
    //which stores for each pattern all possible barcodes, number of mismatches etc.
    std::vector<std::pair<std::string, char> > pattern = this->generate_barcode_patterns(input);

    //create output files and write headers for demultiplexed barcodes
    initialize_output_files(input, pattern);

    //run mapping
    this->run_mapping(input);

    //write the barcodes, failed lines, statistics (mismatches per barcode)
    write_file(input, this->get_demultiplexed_reads());
}

template class MappingAroundLinker<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromFastqFilePolicy>;
template class MappingAroundLinker<MapAroundConstantBarcodesAsAnchorPolicy, ExtractLinesFromTxtFilesPolicy>;
