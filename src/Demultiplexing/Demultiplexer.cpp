#include "Demultiplexer.hpp"

/**
* @brief function wrapping the demultiplex_read function of the Mapping class to increment and decrement a counter of elements
* in the current queue. Used to allow processing of only a few lines a time instead of writing all into RAM.
* This function is called from every thread.
**/
/*
template <typename MappingPolicy, typename FilePolicy>
void Demultiplexer<MappingPolicy, FilePolicy>::demultiplex_wrapper(const std::pair<fastqLine, fastqLine>& line,
                                                                    const input& input,
                                                                    std::atomic<int>& lineCount,
                                                                    const unsigned long long& totalReadCount,
                                                                    std::atomic<long long int>& elementsInQueue)
{

    //FOR EVERY BARCODE-PATTERN (GET PATTERNID)
    //try to map read to this pattern until it matches and exit
    bool result = false;

    std::string foundPatternName;
    DemultiplexedLine finalDemultiplexedLine;
    OneLineDemultiplexingStatsPtr finalLineStatsPtr; //result for a single line
    int bestPatternScore = std::numeric_limits<int>::max();

    // map every pattern and save the overall score per pattern 
    // *this->get_barcode_pattern() for global pattern that is shared
    // *(thread_pattern[boost::this_thread::get_id()])

    for(BarcodePatternPtr pattern : *(thread_pattern[boost::this_thread::get_id()]) )
    {
        //score for this specific pattern
        int tmpPatternScore = std::numeric_limits<int>::max();

        //if we save stats initialize the temporary ones here
        OneLineDemultiplexingStatsPtr lineStatsPtr; //result for a single line
        if(input.writeStats){lineStatsPtr = std::make_shared<OneLineDemultiplexingStats>();}
        else{lineStatsPtr = nullptr;}

        //create result object in which we safe the result
        DemultiplexedLine tmpDemultiplexedLine;

        //pattern contains DNA also store the read name
        if(pattern->containsDNA)
        {
            //use ONLY the forward read name (in the DNA/barcode file later we also add threadID and a readID within thread for unique names)
            tmpDemultiplexedLine.readName = line.first.name;
        }

        //write demultiplexed information into demultiplexedLine, this is passed by reference and can be accessed here
        if(this->demultiplex_read(line, tmpDemultiplexedLine, pattern, input, lineCount, totalReadCount, tmpPatternScore, lineStatsPtr) && tmpPatternScore < bestPatternScore)
        {
            //as soon as a pattern matches, we exit and safe it!
            foundPatternName = pattern->patternName;
            result = true;
            finalDemultiplexedLine = tmpDemultiplexedLine;
            finalLineStatsPtr = lineStatsPtr;
            bestPatternScore = tmpPatternScore;
            //TODO: potentailly break here if we do not want to map all patterns
            //break;
        }

        //in case we have only ONE PATTERN we can get statistics for the failed line, otherwise not
        if(this->get_barcode_pattern()->size() == 1)
        {
            finalLineStatsPtr = lineStatsPtr;
        }
    }

    //BARCODE only information is stored (e.g., protein+barcode, guide+barcode)
    if(result && !finalDemultiplexedLine.containsDNA)
    {
        //store in a shared object
        this->fileWriter->add_demultiplexed_line(foundPatternName, finalDemultiplexedLine.barcodeList);
    }
    else if(result && finalDemultiplexedLine.containsDNA)
    {
        //write out immediately into file for thread (bcs. RNA reads are not very repretitive and might take quite some memory)
        // call write_dna_line
        this->fileWriter->write_dna_line(fileWriter->get_streams_for_threadID(boost::this_thread::get_id(), foundPatternName), finalDemultiplexedLine, boost::this_thread::get_id());
    }
    else if(!result && input.writeFailedLines)
    {
        //write failed line to file, get the filestreams for thread
        std::pair<std::shared_ptr<std::ofstream>, std::shared_ptr<std::ofstream>> failedFileStream = fileWriter->get_failedStream_for_threadID_at(boost::this_thread::get_id());  
        fileWriter->write_failed_line(failedFileStream, line);
    }

    //update statistics
    if(input.writeStats)
    {
        //if we mapped the line store mapping information
        fileWriter->update_stats(finalLineStatsPtr, result, foundPatternName, finalDemultiplexedLine.barcodeList);
    }

    //count down elements to process
    --elementsInQueue;

}*/

template <typename MappingPolicy, typename FilePolicy>
void Demultiplexer<MappingPolicy, FilePolicy>::demultiplex_wrapper_batch(const std::vector<std::pair<fastqLine, fastqLine>>& line_vector,
                                                                    const input& input,
                                                                    std::atomic<unsigned long long>& lineCount,
                                                                    const unsigned long long& totalReadCount,
                                                                    std::atomic<long long int>& elementsInQueue)
{
    //FOR EVERY BARCODE-PATTERN (GET PATTERNID)
    //try to map read to this pattern until it matches and exit
    //this->get_barcode_pattern();  or //(thread_pattern[boost::this_thread::get_id()]) ;
    std::shared_ptr<std::vector<BarcodePatternPtr>> patternVectorPtr = this->get_barcode_pattern();

    //temporary result objects
    //mapping pattern name -> DemultipelxedReads (batch of lines)
    std::unordered_map<std::string, DemultiplexedReads> demultiplexedBatch;
    //initialize result for all patterns
    for(const BarcodePatternPtr& pattern : *patternVectorPtr)
    {
        demultiplexedBatch.emplace(pattern->patternName, DemultiplexedReads(line_vector.size()));
    }
    std::vector<std::pair<fastqLine, fastqLine>> failedLinesBatch;

    for(std::pair<fastqLine, fastqLine> line : line_vector)
    {
        ++lineCount;
        bool result = false;

        std::string foundPatternName;
        DemultiplexedLine finalDemultiplexedLine;
        OneLineDemultiplexingStatsPtr finalLineStatsPtr; //result for a single line
        int bestPatternScore = std::numeric_limits<int>::max();

        //map every pattern and save the overall score per pattern 
        //*this->get_barcode_pattern() for global pattern that is shared
        //*(thread_pattern[boost::this_thread::get_id()])

        for(BarcodePatternPtr pattern : *patternVectorPtr)
        {
            //score for this specific pattern
            int tmpPatternScore = std::numeric_limits<int>::max();

            //if we save stats initialize the temporary ones here
            OneLineDemultiplexingStatsPtr lineStatsPtr; //result for a single line
            if(input.writeStats){lineStatsPtr = std::make_shared<OneLineDemultiplexingStats>();}
            else{lineStatsPtr = nullptr;}

            //create result object in which we safe the result
            DemultiplexedLine tmpDemultiplexedLine;
            //safe the read name - FROM FORWARD READ (later we also add threadID and a readID within thread for unique names)
            tmpDemultiplexedLine.readName = line.first.name;
            
            //write demultiplexed information into demultiplexedLine, this is passed by reference and can be accessed here
            if(this->demultiplex_read(line, tmpDemultiplexedLine, pattern, input, lineCount, totalReadCount, tmpPatternScore, lineStatsPtr) && tmpPatternScore < bestPatternScore)
            {
                //as soon as a pattern matches, we exit and safe it!
                foundPatternName = pattern->patternName;
                result = true;
                finalDemultiplexedLine = tmpDemultiplexedLine;
                finalLineStatsPtr = lineStatsPtr;
                bestPatternScore = tmpPatternScore;
                //TODO: potentailly break here if we do not want to map all patterns
                //break;
            }

            //in case we have only ONE PATTERN we can get statistics for the failed line, otherwise not
            if(this->get_barcode_pattern()->size() == 1)
            {
                finalLineStatsPtr = lineStatsPtr;
            }
        }

        //push demultiplexed line into vector for the right pattern, after processing batch it will be written
        if(result)
        {
            demultiplexedBatch.at(foundPatternName).store_demultiplexed_read(finalDemultiplexedLine);
        }
        else if(!result && input.writeFailedLines)
        {
            failedLinesBatch.emplace_back(line);
        }

        if(input.writeStats)
        {
            //store mapping information
            std::shared_ptr<DemultiplexingStats> threadTmpStats = fileWriter->get_demultiplexing_stats(boost::this_thread::get_id());
            fileWriter->update_stats(threadTmpStats, finalLineStatsPtr, result, foundPatternName, finalDemultiplexedLine.barcodeList);
        }

        //BARCODE only information is stored (e.g., protein+barcode, guide+barcode)
        /*if(result && !finalDemultiplexedLine.containsDNA)
        {
            //store in a shared object
            this->fileWriter->add_demultiplexed_line(foundPatternName, finalDemultiplexedLine.barcodeList);
        }
        else if(result && finalDemultiplexedLine.containsDNA)
        {
            //write out immediately into file for thread (bcs. RNA reads are not very repretitive and might take quite some memory)
            // call write_dna_line
            this->fileWriter->write_dna_line(fileWriter->get_streams_for_threadID(boost::this_thread::get_id(), foundPatternName), finalDemultiplexedLine, boost::this_thread::get_id());
        }
        else if(!result && input.writeFailedLines)
        {
            //write failed line to file, get the filestreams for thread
            std::pair<std::shared_ptr<std::ofstream>, std::shared_ptr<std::ofstream>> failedFileStream = fileWriter->get_failedStream_for_threadID_at(boost::this_thread::get_id());  
            fileWriter->write_failed_line(failedFileStream, line);
        }*/

        //update statistics
       // if(input.writeStats)
        //{
            //if we mapped the line store mapping information
        //    fileWriter->update_stats(finalLineStatsPtr, result, foundPatternName, finalDemultiplexedLine.barcodeList);
       // }

        //count down elements to process
        --elementsInQueue;
    }   

    //NEW BATCH WRITING FUNCTIONS
    //write found lines
    for(const BarcodePatternPtr& pattern : *patternVectorPtr)
    {
        fileWriter->write_demultiplexed_batch(fileWriter->get_streams_for_threadID(boost::this_thread::get_id(), pattern->patternName), 
                                              demultiplexedBatch.at( pattern->patternName), 
                                              boost::this_thread::get_id(), pattern->containsDNA);  
    }
    
    //write failed lines   
    if(input.writeFailedLines)
    {
        std::pair<std::shared_ptr<std::ofstream>, std::shared_ptr<std::ofstream>> failedFileStream = fileWriter->get_failedStream_for_threadID_at(boost::this_thread::get_id());  
        fileWriter->write_failed_lines(failedFileStream, failedLinesBatch);
    }

}

/// overwritten run_mapping function to allow processing of only a subset of fastq lines at a time
template <typename MappingPolicy, typename FilePolicy>
void Demultiplexer<MappingPolicy, FilePolicy>::run_mapping(const input& input)
{
    //generate a pool of threads
    boost::asio::thread_pool pool(input.threads); //create thread pool
    //initialize thread-dependent tmp files
    fileWriter->initialize_thread_streams(pool, input.threads);
    //initialize the statistics
    if(input.writeStats)
    {
        //statitstics need the barcode pattern to initialize, e.g., MM per pattern at each position
        fileWriter->initialize_statistiscs(this->get_barcode_pattern());
    }
    //initialize copies of Mapping pattern for all threads - reading patterns seems to be no bottleneck
    //share pattern object across threads
    //initialize_thread_patterns(pool, input.threads);

    //read line by line and add to thread pool
    this->FilePolicy::init_file(input.inFile, input.reverseFile);
    std::pair<fastqLine, fastqLine> line;
    std::atomic<unsigned long long> lineCount = 0; //using atomic<int> as thread safe read count
    std::atomic<long long int> elementsInQueue(0);
    unsigned long long totalReadCount = FilePolicy::get_read_number();

    /*while(FilePolicy::get_next_line(line))
    {
        //wait to enqueue new elements in case we have a maximum bucket size
        if(input.fastqReadBucketSize>0)
        {
            while(input.fastqReadBucketSize <= elementsInQueue.load()){}
        }
        //increase job count and push the job in the queue
        
        ++elementsInQueue;
        boost::asio::post(pool, std::bind(&Demultiplexer::demultiplex_wrapper, this, line, input, ++lineCount, totalReadCount, std::ref(elementsInQueue)));
    }*/

    std::cout << "START DEMULTIPLEXING\n";
    //btch processing of lines (better for thread scheduling)
    // by default we parse 100K lines, making a batch size of 100 lines for 10 threads
    // benchmarked for a factor of 1, 10, 100. 100 gave best performance
    unsigned long batchSize = input.fastqReadBucketSize/ (input.threads * 100);
    std::vector<std::pair<fastqLine, fastqLine>> lineBatch;
    lineBatch.reserve(batchSize); // e.g., batchSize = 100
    while (FilePolicy::get_next_line(line)) 
    {
        lineBatch.push_back(std::move(line));
        ++elementsInQueue;

        if (lineBatch.size() == batchSize) 
        {
            // wait if we exceed the queue size limit
            if (input.fastqReadBucketSize > 0) 
            {
                while (input.fastqReadBucketSize <= elementsInQueue.load()) {}
            }

            boost::asio::post(pool, std::bind(
                &Demultiplexer::demultiplex_wrapper_batch,
                this,
                std::move(lineBatch),
                input,
                std::ref(lineCount), // starting line number
                totalReadCount,
                std::ref(elementsInQueue)
            ));

            lineBatch.clear();
            lineBatch.reserve(batchSize);
        }
    }
    // Handle any remaining lines
    if (!lineBatch.empty()) {
        ++elementsInQueue;
        boost::asio::post(pool, std::bind(
            &Demultiplexer::demultiplex_wrapper_batch,
            this,
            std::move(lineBatch),
            input,
            std::ref(lineCount),
            totalReadCount,
            std::ref(elementsInQueue)
        ));
    }

    //iterate through the barcode map and let threads 
    pool.join();

    printProgress(1); std::cout << "\n"; // end the progress bar
    if(totalReadCount != ULLONG_MAX && input.writeStats)
    {
        //combine statistics before writing results
        std::vector<std::shared_ptr<DemultiplexingStats>> statsList;
        for (const std::pair<boost::thread::id, std::shared_ptr<DemultiplexingStats>> threadStatPair : this->fileWriter->get_statistics_map()) 
        {
            statsList.push_back(threadStatPair.second);
        }
        this->fileWriter->combine_statistics(statsList);
        std::cout << "=>\tPERFECT MATCHES: " << std::to_string((unsigned long long)(100*(this->fileWriter->get_perfect_matches())/(double)totalReadCount)) 
                  << "% | MODERATE MATCHES: " << std::to_string((unsigned long long)(100*(this->fileWriter->get_moderat_matches())/(double)totalReadCount))
                  << "% | MISMATCHES: " << std::to_string((unsigned long long)(100*(this->fileWriter->get_failed_matches())/(double)totalReadCount)) << "%\n";
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
    std::cout << "Initializing barcode pattern\n";
    this->generate_barcode_patterns(input);

    //create output files and write headers for demultiplexed barcodes
    fileWriter = std::make_shared<DemultiplexedResult>(DemultiplexedResult(input, this->get_barcode_pattern()));

    //create empty dict for mismatches per barcode
    //THIS SHould noW GET iniTIalized in OUTPUTFILEWRITER in initialieStats
    //if(input.writeStats)
    //{
    //    this->initializeStats();
    //}

    //run mapping
    this->run_mapping(input);

    //write the barcodes, failed lines, statistics (mismatches per barcode)
    //TODO:
    //iterate over the different result maps for the various barcodePatterns

    //write all final files
    //instead write now batches periodically
    fileWriter->write_output(input);
}

template class Demultiplexer<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy>;
template class Demultiplexer<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy>;
template class Demultiplexer<MapEachBarcodeSequentiallyPolicyPairwise, ExtractLinesFromFastqFilePolicyPairedEnd>;