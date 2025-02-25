#include "Demultiplexer.hpp"

//initialize the statistics file/ lines that could not be mapped
//FILES: mismatches per barcode / mismatches per barcodePattern/ failedLines
void OutputFileWriter::initialize_additional_output(std::string output)
{
    //remove output
    barcodeMismatches = output + "/BarcodeMismatches.txt";
    patternMismatches = output + "/PatternMismatches.txt";
    failedLines = output + "/FailedLines.txt";

    // remove outputfile if it exists
    std::remove(barcodeMismatches.c_str());
    std::remove(patternMismatches.c_str());
    std::remove(failedLines.c_str());

    std::ofstream barcodeMMFile(barcodeMismatches.c_str());
    std::ofstream patternMMFile(patternMismatches.c_str());
    std::ofstream filedLinesFile(failedLines.c_str());
}

//this function only initialized the output files that are needed for a specific pattern
//like fastq, demultiplexed-read files
void OutputFileWriter::initialize_output_for_pattern(std::string output, const BarcodePatternPtr pattern)
{
    // 1.) INITIALIZE TWO FILES
    FinalPatternFiles patternOutputs;
    //txt-file with barcodes
    patternOutputs.barcodeFile = output + "/" + pattern->patternName + ".txt";
    std::remove(patternOutputs.barcodeFile.c_str());
    //fastq-file with the RNA sequence
    patternOutputs.dnaFile = "";
    if(pattern->containsDNA)
    {
        patternOutputs.dnaFile = output + "/" + pattern->patternName + ".fastq";
        std::remove(patternOutputs.dnaFile.c_str());
    }

    // 2.) STORE OUTPUT-FILE NAMES IN LIST
    finalFiles[pattern->patternName] = patternOutputs;

    // 3.) WRITE POTENTIAL HEADER OF BARCODE FILE
    //write header line for barcode file: e.g.: [ACGGCATG][BC1.txt][15X]
    std::ofstream barcodeOutputStream;   
    barcodeOutputStream.open(patternOutputs.barcodeFile, std::ofstream::app);

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

/// calls output initializer functions and gets the barcode mapping structure from Mapping object, since this will the header of the output file
// create backbone files for barcoding patterns that will be mapped: e.g.: FASTQ for RNA, txt with heads for barcode-files for CI, spatial, other stuff
//the file will be anmed after pattern name
void OutputFileWriter::initialize_output_files(const input& input, const MultipleBarcodePatternVectorPtr& barcodePatternList)
{
    //TO DO
    //parse through the barcodePatterns, make file of pattern name
    for(const BarcodePatternPtr barcodePattern : *barcodePatternList)
    {
        initialize_output_for_pattern(input.outPath, barcodePattern);
    }

    //create universal output files
    // 2 statistics files: mismatches per pattern, and mismatches in barcodes
    // file with failed lines
    initialize_additional_output(input.outPath);
}

void OutputFileWriter::initialize_tmp_file()
{

    //create a map of pattern name to open streams for barcodes, fastq for every thread
    //those thread files have no header, the header is only written in final file
    std::unordered_map<std::string, TmpPatternStream> tmpFileStreams;
    //the files that r written immediately are failed/ DNA files
    for(auto fileIt = finalFiles.begin(); fileIt != finalFiles.end(); ++fileIt)
    {
        TmpPatternStream tmpStream;
        // Create an ofstream pointer and open the file
        std::ofstream* outFileBarcode = new std::ofstream(fileIt->second.barcodeFile);
        if (!outFileBarcode->is_open()) 
        {
            std::cerr << "Error opening file: " << fileIt->second.barcodeFile << std::endl;
            exit(EXIT_FAILURE);
        }
        tmpStream.barcodeStream = outFileBarcode;
        if(fileIt->second.dnaFile != "")
        {
            std::ofstream* outFileDna = new std::ofstream(fileIt->second.dnaFile);
            if (!outFileDna->is_open()) 
            {
                std::cerr << "Error opening file: " << fileIt->second.dnaFile << std::endl;
                exit(EXIT_FAILURE);
            }
            tmpStream.dnaStream = outFileDna;
        }
        tmpFileStreams[fileIt->first] = tmpStream;
    }

    std::unique_lock<std::mutex> fileLock(*threadFileOpenerMutex);
    --(*threadToInitializePtr);
    //add the new list of streams to map
    tmpStreamMap[boost::this_thread::get_id()] = tmpFileStreams;
    //decrease number of threads that need initialization, when all are initialized we can continue program in main function
    if (*threadToInitializePtr == 0) 
    {
        cvPtr->notify_all();  // Notify when all tasks are done
        fileLock.unlock();
    }
    else
    {
        //lock threads if we did not finish them - like this we can be 100% sure that evey thread gets called EXACTLY ONCE
        // for initialization
        fileLock.unlock();
        //wait until all threads are here
        std::unique_lock<std::mutex> lock(*threadWaitingMutex);
        cvPtr->wait(lock, [&] { return (*threadToInitializePtr) == 0; });
        lock.unlock();
    }

}

void OutputFileWriter::initialize_thread_streams(boost::asio::thread_pool& pool, const int threadNum)
{
    //add the initialiozation fucntion thread times, where every thread waits until threadsStarted is equal
    //the number of threads to initialize every thread exactly once
    for(int i = 0; i < threadNum; ++i)
    {
        boost::asio::post(pool, std::bind(&OutputFileWriter::initialize_tmp_file, this));
    }

    //only continue once all htreads have finished 
    // (we can NOT join, to not change the thread_ids, and waiting within thread only might let threadToInitialize go out of scope)
    std::unique_lock<std::mutex> lock(*threadWaitingMutex);
    cvPtr->wait(lock, [&] { return (*threadToInitializePtr) == 0; });
    lock.unlock();

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










/**
* @brief function wrapping the demultiplex_read function of the Mapping class to increment and decrement a counter of elements
* in the current queue. Used to allow processing of only a few lines a time instead of writing all into RAM.
**/
template <typename MappingPolicy, typename FilePolicy>
void Demultiplexer<MappingPolicy, FilePolicy>::demultiplex_wrapper(std::pair<const std::string&, const std::string&> line,
                                                            const input& input,
                                                            const unsigned long long lineCount,
                                                            const unsigned long long& totalReadCount,
                                                            std::atomic<long long int>& elementsInQueue)
{

    //create result object in which we safe the result
    DemultiplexedLine demultiplexedLine;

    //FOR EVERY BARCODE-PATTERN (GET PATTERNID)
    bool result = false;
    for(BarcodePatternPtr pattern : *this->get_barcode_pattern())
    {
        if(pattern->containsDNA)
        {
            // set the quality of the read and also the read name

        }

        //write demultiplexed information into demultiplexedLine, this is passed by reference and can be accessed here
        if(this->demultiplex_read(line, demultiplexedLine, pattern, input, lineCount, totalReadCount, false))
        {
            //as soon as a pattern matches, we exit and safe it!
            //we DO NOT check if other patterns match as well -> responsibility of user to design non-ambiguous patterns
            break;
        }
    }

    //decide weather demultiplexedLine is stored or written to an output file
    if(result && !demultiplexedLine.containsDNA)
    {
        //store in a shared object

    }
    else if(result && demultiplexedLine.containsDNA)
    {
        //write out immediately into file for thread
        std::cout << boost::this_thread::get_id() << "\n";

    }
    else if(!result && input.writeFailedLines)
    {
        //write failed line to file
        //write_failed_line(input, line);
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
    std::pair<std::string, std::string> line;
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

    fileWriter->write_output();

}

template class Demultiplexer<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy>;
template class Demultiplexer<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy>;
template class Demultiplexer<MapEachBarcodeSequentiallyPolicyPairwise, ExtractLinesFromFastqFilePolicyPairedEnd>;