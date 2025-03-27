#include "Demultiplexer.hpp"

/**
* @brief function wrapping the demultiplex_read function of the Mapping class to increment and decrement a counter of elements
* in the current queue. Used to allow processing of only a few lines a time instead of writing all into RAM.
**/
template <typename MappingPolicy, typename FilePolicy>
void Demultiplexer<MappingPolicy, FilePolicy>::demultiplex_wrapper(const std::pair<fastqLine, fastqLine>& line,
                                                            const input& input,
                                                            const unsigned long long lineCount,
                                                            const unsigned long long& totalReadCount,
                                                            std::atomic<long long int>& elementsInQueue)
{

    //create result object in which we safe the result
    DemultiplexedLine demultiplexedLine;

    //FOR EVERY BARCODE-PATTERN (GET PATTERNID)
    bool result = false;
    std::string foundPatternName;
    for(BarcodePatternPtr pattern : *this->get_barcode_pattern())
    {
        //pattern contains DNA also store the read name
        if(pattern->containsDNA)
        {
            //use ONLY the forward read name (in the DNA/barcode file later we also add threadID and a readID within thread for unique names)
            demultiplexedLine.dnaName = line.first.name;
        }

        //write demultiplexed information into demultiplexedLine, this is passed by reference and can be accessed here
        if(this->demultiplex_read(line, demultiplexedLine, pattern, input, lineCount, totalReadCount, false))
        {
            //as soon as a pattern matches, we exit and safe it!
            //we DO NOT check if other patterns match as well -> responsibility of user to design non-ambiguous patterns
            foundPatternName = pattern->patternName;
            result = true;
            break;
        }
    }

    //BARCODE only information is stored (e.g., protein+barcode, guide+barcode)
    if(result && !demultiplexedLine.containsDNA)
    {
        //store in a shared object
        this->fileWriter->add_demultiplexed_line(foundPatternName, demultiplexedLine.barcodeList);
    }
    else if(result && demultiplexedLine.containsDNA)
    {
        //write out immediately into file for thread (bcs. RNA reads are not very repretitive and might take quite some memory)
        // call write_dna_line
        this->fileWriter->write_dna_line(fileWriter->get_streams_for_threadID(boost::this_thread::get_id(), foundPatternName), demultiplexedLine, boost::this_thread::get_id());
    }
    else if(!result && input.writeFailedLines)
    {
        //write failed line to file, get the filestreams for thread
        std::pair<std::shared_ptr<std::ofstream>, std::shared_ptr<std::ofstream>> failedFileStream = fileWriter->get_failedStream_for_threadID_at(boost::this_thread::get_id());
        fileWriter->write_failed_line(failedFileStream, line);

    }
    --elementsInQueue;

}

/// overwritten run_mapping function to allow processing of only a subset of fastq lines at a time
template <typename MappingPolicy, typename FilePolicy>
void Demultiplexer<MappingPolicy, FilePolicy>::run_mapping(const input& input)
{
    std::cout << "START DEMULTIPLEXING\n";

    //generate a pool of threads
    boost::asio::thread_pool pool(input.threads); //create thread pool
    //initialize thread-dependent tmp files
    fileWriter->initialize_thread_streams(pool, input.threads);

    //read line by line and add to thread pool
    this->FilePolicy::init_file(input.inFile, input.reverseFile);
    std::pair<fastqLine, fastqLine> line;
    unsigned long long lineCount = 0; //using atomic<int> as thread safe read count
    std::atomic<long long int> elementsInQueue(0);
    unsigned long long totalReadCount = FilePolicy::get_read_number();

    while(FilePolicy::get_next_line(line))
    {
        //wait to enqueue new elements in case we have a maximum bucket size
        if(input.fastqReadBucketSize>0)
        {
            while(input.fastqReadBucketSize <= elementsInQueue.load()){}
        }
        //increase job count and push the job in the queue
        ++elementsInQueue;
        boost::asio::post(pool, std::bind(&Demultiplexer::demultiplex_wrapper, this, line, input, ++lineCount, totalReadCount, std::ref(elementsInQueue)));
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
void Demultiplexer<MappingPolicy, FilePolicy>::run(const input& input)
{

    //from the basic information within patterns generate a more complex barcodePattern object
    //which stores for each pattern all possible barcodes, number of mismatches etc.
    this->generate_barcode_patterns(input);

    //create output files and write headers for demultiplexed barcodes
    fileWriter = std::make_shared<OutputFileWriter>(OutputFileWriter(input, this->get_barcode_pattern()));

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

    //write all final files
    fileWriter->write_output(input);
}

template class Demultiplexer<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy>;
template class Demultiplexer<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy>;
template class Demultiplexer<MapEachBarcodeSequentiallyPolicyPairwise, ExtractLinesFromFastqFilePolicyPairedEnd>;