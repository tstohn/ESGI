#include "DemultiplexedResult.hpp"

void DemultiplexedResult::concatenateFiles(const std::vector<std::string>& tmpFileList, const std::string& outputFile) 
{
    for (const std::string& file : tmpFileList) 
    {
        std::ifstream in(file, std::ios::in| std::ios::binary);  // Read written files
        if (!in) 
        {
            std::cerr << "Error opening input file: " << file << "\n";
            continue;
        }

        //get content as string
        std::ostringstream buffer;
        buffer << in.rdbuf();  
        std::string content = buffer.str();
        in.close();

        //write temporary file into final file
        std::ofstream out(outputFile, std::ios::app | std::ios::binary);  // Final output file
        out.write(content.c_str(), content.size());
        out.close();

        // Delete the temporary file
        if (std::remove(file.c_str()) != 0) 
        {
            std::cerr << "Error: Could not delete " << file << std::endl;
        }
    }
}

void DemultiplexedResult::update_stats(std::shared_ptr<DemultiplexingStats> threadTmpStats, OneLineDemultiplexingStatsPtr lineStatsPtr, 
                                       bool result, std::string& foundPatternName, std::vector<std::string>& barcodeList)
{
    threadTmpStats->update(lineStatsPtr, result, foundPatternName, barcodeList);
}

//writes the dna (fastq) and barcode (tsv) data
void DemultiplexedResult::write_dna_line(TmpPatternStream& dnaLineStream, std::ostringstream& lineBuffer, const DemultiplexedLine& demultiplexedLine, const boost::thread::id& threadID)
{
    //write RNA data to dnaBuffer (FASTQ)
    lineBuffer << "@" << threadID << "_" << ++dnaLineStream.lineNumber << "_" << demultiplexedLine.readName << "\n";

    lineBuffer << demultiplexedLine.dna << "\n";
    lineBuffer << "+\n";
    lineBuffer << demultiplexedLine.dnaQuality << "\n";
}

void DemultiplexedResult::write_barcode_line(TmpPatternStream& barcodeLineStream, std::ostringstream& lineBuffer, 
                                             const DemultiplexedLine& demultiplexedLine, const boost::thread::id& threadID)
{
    // Build the line in a temporary buffer
    lineBuffer << threadID << "_" << ++barcodeLineStream.lineNumber << "_" << demultiplexedLine.readName << "\t";
    for (size_t i = 0; i < demultiplexedLine.barcodeList.size(); ++i) 
    {
        lineBuffer << demultiplexedLine.barcodeList[i];
        if (i + 1 < demultiplexedLine.barcodeList.size()) 
        {
            lineBuffer << "\t";
        } else {
            lineBuffer << "\n";
        }
    }
}

void DemultiplexedResult::write_demultiplexed_batch(TmpPatternStream& lineStream, const DemultiplexedReads& demultiplexedBatch, const boost::thread::id& threadID, const bool containsDNA)
{
    const std::vector<DemultiplexedLine>& lines = demultiplexedBatch.get_all_reads();

    //write barcodes to temporary files
    std::ostringstream barcodeLineBuffer; //temporary buffer for the batch of demultiplexed lines
    lineStream.lineNumber = 0; //reset line number, it is incremented while writing lines
    for(const DemultiplexedLine& line : lines)
    {
        //local buffer to save the whole batch of demultiplexed lines
        write_barcode_line(lineStream, barcodeLineBuffer, line, threadID);
    }
    // Write the full line in one go (one I/O call)
    *(lineStream.barcodeStream) << barcodeLineBuffer.str();

    //if dna is present write it to seperate file
    if(containsDNA)
    {
        std::ostringstream dnaLineBuffer;
        lineStream.lineNumber = 0; //reset line number, it is incremented while writing lines
        for(const DemultiplexedLine& line : lines)
        {
            write_dna_line(lineStream, dnaLineBuffer, line, threadID);
        }
        *(lineStream.dnaStream) << dnaLineBuffer.str();
    }
}

void DemultiplexedResult::close_and_concatenate_fileStreams(const input& input)
{

    //1.) CLOSE all file streams
    //iterate over all threads
    for (auto threadID_patternStreamMap_pair : tmpStreamMap) 
    {
        //CLOSED FAILED LINE STREAMS: failed file per thread
        std::pair<std::shared_ptr<std::ofstream>, std::shared_ptr<std::ofstream>> failedFileStream = get_failedStream_for_threadID_at(threadID_patternStreamMap_pair.first);
        //close the first stream definitely
        failedFileStream.first->close();
        if(failedFileStream.second != nullptr)
        {
            //and the second only if we have paired-end sequencing
            failedFileStream.second->close();
        }

        //CLOSE DNA-LINE STREAMS: for failed (DNA, barcode-tsv)-pair we have one file per pattern (with DNA)
        for (auto patternStreams_pair : threadID_patternStreamMap_pair.second) 
        {
            //close barcode stream if existing
            std::shared_ptr<std::ofstream> barcodeStream = get_barcodeStream_for_threadID_at(threadID_patternStreamMap_pair.first, patternStreams_pair.first);
            barcodeStream->close();

            //close DNA stream if existing
            std::shared_ptr<std::ofstream> dnaStream = get_dnaStream_for_threadID_at(threadID_patternStreamMap_pair.first, patternStreams_pair.first);
            if(dnaStream != nullptr)
            {
                dnaStream->close();
            }
        }
    }

    //2.) CONCATENATE all files
    //CONCATENATE FAILED LINES
    //FAILED LINES FORWARDS
    std::vector<std::string> failedFileListFW;
    //ALWAYS concatenate the failed lines for the first fileFilesName (FW or single-read)
    size_t dotPos = failedLines.first.find_last_of('.');  // Find the last dot
    std::string failedLinesTmpFileNameFW;
    for(int i = 0; i < input.threads; ++i)
    {
        failedLinesTmpFileNameFW = failedLines.first.substr(0, dotPos) + std::to_string(i) + failedLines.first.substr(dotPos);
        failedFileListFW.push_back(failedLinesTmpFileNameFW);
    }
    concatenateFiles(failedFileListFW, failedLines.first);
    //FAILED-LINES REVERSE (if we have paired-end also concatenated RV reads)
    if(failedLines.second != "")
    {
        std::vector<std::string> failedFileListRV;
        dotPos = failedLines.second.find_last_of('.');  // Find the last dot
        std::string failedLinesTmpFileNameRV;
        for(int i = 0; i < input.threads; ++i)
        {
            failedLinesTmpFileNameRV = failedLines.second.substr(0, dotPos) + std::to_string(i) + failedLines.second.substr(dotPos);
            failedFileListRV.push_back(failedLinesTmpFileNameRV);
        }
        concatenateFiles(failedFileListRV, failedLines.second);
    }

    //CONCATENATED DEMULTIPLEXED-LINES-FILES
    for(auto fileIt = finalFiles.begin(); fileIt != finalFiles.end(); ++fileIt)
    {
        //only create a temporary stream for a pattern that does contain a DNA region
        if(fileIt->second.dnaFile != "")
        {
            //temporary file vectors
            std::vector<std::string> dnaFileList; // fastq files of,e.g., RNA/DNA
            std::vector<std::string> barcodeFileList; //tsv files for corersponding fastq files
        
            //list all temporary dna-files (FASTQ)/ barcode-files (tsv) and combine them
            for(int i = 0; i < input.threads; ++i)
            {
                //pattern and thread specific DNA file
                size_t dotPos = fileIt->second.dnaFile.find_last_of('.');  // Find the last dot
                std::string dnaTmpFileName = fileIt->second.dnaFile.substr(0, dotPos) + std::to_string(i) + fileIt->second.dnaFile.substr(dotPos);
                dnaFileList.push_back(dnaTmpFileName);

                //pattern and thread specific barcode file
                dotPos = fileIt->second.barcodeFile.find_last_of('.');  // Find the last dot
                std::string barcodeTmpFileName = fileIt->second.barcodeFile.substr(0, dotPos) + std::to_string(i) + fileIt->second.barcodeFile.substr(dotPos);
                barcodeFileList.push_back(barcodeTmpFileName);
            }

            concatenateFiles(dnaFileList, fileIt->second.dnaFile);
            concatenateFiles(barcodeFileList, fileIt->second.barcodeFile);
        }
        else //concatenate barcode only files
        {
            //temporary file vectors
            std::vector<std::string> barcodeFileList; //tsv files for corersponding fastq files
        
            for(int i = 0; i < input.threads; ++i)
            {
                //pattern and thread specific barcode file
                dotPos = fileIt->second.barcodeFile.find_last_of('.');  // Find the last dot
                std::string barcodeTmpFileName = fileIt->second.barcodeFile.substr(0, dotPos) + std::to_string(i) + fileIt->second.barcodeFile.substr(dotPos);
                barcodeFileList.push_back(barcodeTmpFileName);
            }
            concatenateFiles(barcodeFileList, fileIt->second.barcodeFile);
        }
    }
}


void DemultiplexedResult::write_output(const input& input)
{
    //if tmp files were written (for failed/ DNA&barcode reads)
    close_and_concatenate_fileStreams(input);

    //write statistics output (already concatenated when demultiplexing ends in demultiplexer)
    if(input.writeStats)
    {
        //write combined statistics (finalFiles.size() just tells the function if we have one or more patterns, to write until
        //where we could map for failed lines)
        finalStat.write(input.outPath, input.prefix, finalFiles.size());
    }

    //write all barcode-only files (data is still in memory at this point)
    //NO LONGER NEEDED to write at end: we concatenate the tmp-files of threads
    /*for (const auto& [patternName, demultiplexedReadsPtrTmp] : demultiplexedReads) 
    {
        write_demultiplexed_barcodes(input, demultiplexedReadsPtrTmp->get_all_reads(), stripQuotes(patternName));
    }*/
    
    //write the results
    /*if(input.writeStats)
    {
        dxStat.write(input.outPath, input.prefix, demultiplexedReads.size());
    }*/
}

//initialize the statistics file/ lines that could not be mapped
//FILES: mismatches per barcode / mismatches per barcodePattern/ failedLines
void DemultiplexedResult::initialize_additional_output(const input& input)
{

    //we need to initialize 1 or 2 files for failed lines - depending on paired/ single read
    if(input.reverseFile == "")
    {
        std::string failedLinesName = "FailedLines.txt";
        if(input.prefix != "")
        {
            failedLinesName = input.prefix + "_" + failedLinesName;
        }
        failedLines.first = input.outPath + "/" + failedLinesName;
        failedLines.second = "";

        std::remove(failedLines.first.c_str());
        std::ofstream failedLineFileFW(failedLines.first.c_str(), std::ios::out | std::ios::binary);
        if (!failedLineFileFW) 
        {
            std::cerr << "Error: Could not create " << failedLines.first << "!\n";
        }
        failedLineFileFW.close();  // Close the file
    }
    else
    {
        std::string failedLinesFWName = "FailedLines_FW.txt";
        if(input.prefix != "")
        {
            failedLinesFWName = input.prefix + "_" + failedLinesFWName;
        }
        failedLines.first = input.outPath + "/" + failedLinesFWName;
        std::remove(failedLines.first.c_str());
        std::ofstream failedLineFileFW(failedLines.first.c_str(), std::ios::out | std::ios::binary);
        if (!failedLineFileFW) 
        {
            std::cerr << "Error: Could not create " << failedLines.first << "!\n";
        }
        failedLineFileFW.close();  // Close the file

        std::string failedLinesRVName = "FailedLines_RV.txt";
        if(input.prefix != "")
        {
            failedLinesRVName = input.prefix + "_" + failedLinesRVName;
        }
        failedLines.second = input.outPath + "/" + failedLinesRVName;
        std::remove(failedLines.second.c_str());
        std::ofstream failedLineFileRV(failedLines.second.c_str(), std::ios::out | std::ios::binary);
        if (!failedLineFileRV) 
        {
            std::cerr << "Error: Could not create " << failedLines.second << "!\n";
        }
        failedLineFileRV.close();  // Close the file
        
    }

    // std::ofstream barcodeMMFile(barcodeMismatches.c_str());
    // std::ofstream patternMMFile(patternMismatches.c_str());
}

//this function only initialized the output files that are needed for a specific pattern
//like fastq, demultiplexed-read files
void DemultiplexedResult::initialize_output_for_pattern(const std::string& output, const std::string& prefix, const BarcodePatternPtr pattern)
{
    // 1.) INITIALIZE TWO FILES
    FinalPatternFiles patternOutputs;
   
    //fastq-file with the RNA sequence
    patternOutputs.dnaFile = "";
    //if pattern contains DNA we write a tsv and FASTQ file
    if(pattern->containsDNA)
    {
        //tsv-file with barcodes
        std::string barcodeTsvFileName = stripQuotes(pattern->patternName) + ".tsv";
        if(prefix != "")
        {
            barcodeTsvFileName = prefix + "_" + barcodeTsvFileName;
        }
        patternOutputs.barcodeFile = output + "/" + barcodeTsvFileName;
        std::remove(patternOutputs.barcodeFile.c_str());

        //fastq file
        std::string fastqFileName = stripQuotes(pattern->patternName) + ".fastq";
        if(prefix != "")
        {
            fastqFileName = prefix + "_" + fastqFileName;
        }
        patternOutputs.dnaFile = output + "/" + fastqFileName;
        std::remove(patternOutputs.dnaFile.c_str());

        //STORE OUTPUT-FILE NAMES IN LIST
        finalFiles[pattern->patternName] = patternOutputs;
    }
    else //otherwise we write ONLY a tsv file of barcodes and initialize the DemultiplexedReads structure to store found reads
    {
        std::string barcodeTsvFileName = stripQuotes(pattern->patternName) + ".tsv";
        if(prefix != "")
        {
            barcodeTsvFileName = prefix + "_" + barcodeTsvFileName;
        }
        patternOutputs.barcodeFile = output + "/" + barcodeTsvFileName;
        std::remove(patternOutputs.barcodeFile.c_str());

        //for every pattern create a sharedPtr of DemultiplexedReads
        //demultiplexedReads.emplace(pattern->patternName, std::make_shared<DemultiplexedReads>());
        
        //STORE OUTPUT-FILE NAMES IN LIST
        finalFiles[pattern->patternName] = patternOutputs;
    }

    //2) WRITE HEADER OF BARCODE DATA   (if barcodes correspont to a fastq, the first column stores the read names)
    //write header line for barcode file: e.g.: [ACGGCATG][BC1.txt][15X]
    std::ofstream barcodeOutputStream;   
    barcodeOutputStream.open(patternOutputs.barcodeFile, std::ios::out | std::ios::binary);

    //store the read name in the first column
    barcodeOutputStream << "READNAME";
    //print header for the general pattern
    for(size_t bidx = 0; bidx < (pattern->barcodePattern)->size(); ++bidx)
    {
        BarcodePtr bptr = (pattern->barcodePattern)->at(bidx);
        //stop and DNA pattern should not be written
        if( (bptr->name != "*") && (bptr->name != "DNA") && (bptr->name != "-"))
        {
            
            //get the short name for the barcode (instead of whole path) if it copntains a slash
            std::string filename;
            if (bptr->name.find('/') != std::string::npos || bptr->name.find('\\') != std::string::npos) 
            {
                filename = std::filesystem::path(bptr->name).filename().string();
            } else 
            {
                filename = bptr->name;
            }

            //write tabs before next header-col
            //like this we avoid wrong tabs in the front or back by just checking the for first/last headers since some headers aren t written
            //like *, -, ...
            strip_crlf(filename);
            barcodeOutputStream << "\t" << filename;
        }
    }
    //if we are in independent mode (seperate mapping of reverse and forward read, both 5'->3' direction)
    //the headers are stored in seperate vector
    // IMPORTANT: the [-] pattern-elemnt is ALWAYS ONLY stored in the forward read
    if(pattern->independentReversePattern)
    {
        for(size_t bidx = 0; bidx < (pattern->independentReversePattern)->size(); ++bidx)
        {
            BarcodePtr bptr = (pattern->independentReversePattern)->at(bidx);
            //stop and DNA pattern should not be written

            if( (bptr->name != "*") && (bptr->name != "DNA") && (bptr->name != "-"))
            {

                //get the short name for the barcode (instead of whole path) if it copntains a slash
                std::string filename;
                if (bptr->name.find('/') != std::string::npos || bptr->name.find('\\') != std::string::npos) 
                {
                    filename = std::filesystem::path(bptr->name).filename().string();
                } else 
                {
                    filename = bptr->name;
                }

                barcodeOutputStream << "\t" <<  filename;
            }
        }
    }
    barcodeOutputStream << "\n";
    barcodeOutputStream.close();
}

/// calls output initializer functions and gets the barcode mapping structure from Mapping object, since this will the header of the output file
// create backbone files for barcoding patterns that will be mapped: e.g.: FASTQ for RNA, txt with heads for barcode-files for CI, spatial, other stuff
//the file will be anmed after pattern name
void DemultiplexedResult::initialize(const input& input, const MultipleBarcodePatternVectorPtr& barcodePatternList)
{
    //TO DO
    //parse through the barcodePatterns, make file of pattern name
    for(const BarcodePatternPtr& barcodePattern : *barcodePatternList)
    {
        initialize_output_for_pattern(input.outPath, input.prefix, barcodePattern);
    }

    //create universal output files
    // 2 statistics files: mismatches per pattern, and mismatches in barcodes
    // file with failed lines
    initialize_additional_output(input);
}

//increase buffer sizes for faster writing
void DemultiplexedResult::increase_buffer_for_stream(std::shared_ptr<std::ofstream> stream)
{
    std::unique_lock<std::mutex> fileLock(*threadFileOpenerMutex);
    const size_t bufSize = 1 * ONE_MB;
    bufferVector.emplace_back(std::make_unique<char[]>(bufSize));
    stream->rdbuf()->pubsetbuf(bufferVector.back().get(), bufSize);
    fileLock.unlock();
}


//initializes all TEMPORARY FILES for a single thread ID, those files are counted from 0 to THREADNUM-1
// (i is the index that is also used as tmp name suffix - the actual thradID is not used in these names, only in the maps)
void DemultiplexedResult::initialize_tmp_file(const int i)
{

    //create a map of pattern name to open streams for barcodes, fastq for every thread
    //those thread files have no header, the header is only written in final file
    std::unordered_map<std::string, TmpPatternStream> tmpFileStreams;
    size_t dotPos; // temporary variable storing endpoint of file names (position of dot in filename)

    //the files that r written immediately are failed/ DNA files
    // fileIT (first: patternName, second: barcode,DNA,faield file names)
    for(auto fileIt = finalFiles.begin(); fileIt != finalFiles.end(); ++fileIt)
    {
        // Create an ofstream pointer and open the file
        //tmp-file name is final name + threadID
        //create here a stream for DNA/BARCODE patterns
        if(fileIt->second.dnaFile != "")
        {
            TmpPatternStream tmpStream; // struct storing a stream for barcode(tsv) and DNA(FASTQ)
            //TEMPORARY BARCODE-tsv STREAM
            dotPos = fileIt->second.barcodeFile.find_last_of('.');  // Find the last dot
            std::string barcodeTmpFileName = fileIt->second.barcodeFile.substr(0, dotPos) + std::to_string(i) + fileIt->second.barcodeFile.substr(dotPos);
            std::shared_ptr<std::ofstream> outFileBarcode = std::make_shared<std::ofstream>(barcodeTmpFileName);
            increase_buffer_for_stream(outFileBarcode);
            if (!outFileBarcode->is_open()) 
            {
                std::cerr << "Error opening file: " << fileIt->second.barcodeFile << std::endl;
                exit(EXIT_FAILURE);
            }
            tmpStream.barcodeStream = outFileBarcode;

            //TEMPORARY DNA STREAM
            //tmp-file name is final name + threadID
            dotPos = fileIt->second.dnaFile.find_last_of('.');  // Find the last dot
            std::string dnaTmpFileName = fileIt->second.dnaFile.substr(0, dotPos) + std::to_string(i) + fileIt->second.dnaFile.substr(dotPos);
            std::shared_ptr<std::ofstream> outFileDna = std::make_shared<std::ofstream>(dnaTmpFileName);
            increase_buffer_for_stream(outFileDna);
            if (!outFileDna->is_open()) 
            {
                std::cerr << "Error opening file: " << fileIt->second.dnaFile << std::endl;
                exit(EXIT_FAILURE);
            }
            tmpStream.dnaStream = outFileDna;
            tmpFileStreams[fileIt->first] = tmpStream;
        }
        else //create stream for patterns with barcodes only, dnaStream stays a nullptr
        {
            TmpPatternStream tmpStream; // struct storing a stream for barcode(tsv) and DNA(FASTQ)
            //TEMPORARY BARCODE-tsv STREAM
            dotPos = fileIt->second.barcodeFile.find_last_of('.');  // Find the last dot
            std::string barcodeTmpFileName = fileIt->second.barcodeFile.substr(0, dotPos) + std::to_string(i) + fileIt->second.barcodeFile.substr(dotPos);
            std::shared_ptr<std::ofstream> outFileBarcode = std::make_shared<std::ofstream>(barcodeTmpFileName);
            increase_buffer_for_stream(outFileBarcode);
            if (!outFileBarcode->is_open()) 
            {
                std::cerr << "Error opening file: " << fileIt->second.barcodeFile << std::endl;
                exit(EXIT_FAILURE);
            }
            tmpStream.barcodeStream = outFileBarcode;
            tmpFileStreams[fileIt->first] = tmpStream;
        }
    }

    //TEMPORARY FAILED LINE STREAM
    dotPos = failedLines.first.find_last_of('.');  // Find the last dot
    std::string failedLinesTmpFileName = failedLines.first.substr(0, dotPos) + std::to_string(i) + failedLines.first.substr(dotPos);
    
    std::shared_ptr<std::ofstream> outFileFailedLineFW = std::make_shared<std::ofstream>(failedLinesTmpFileName);
    increase_buffer_for_stream(outFileFailedLineFW);
    if (!outFileFailedLineFW->is_open()) 
    {
        std::cerr << "Error opening file: " << failedLines.first << std::endl;
        exit(EXIT_FAILURE);
    }

    // check if we also have a reverse read that we need to store temporarily
    std::shared_ptr<std::ofstream> outFileFailedLineRV = nullptr;
    if(failedLines.second != "")
    {
        dotPos = failedLines.second.find_last_of('.');  // Find the last dot
        std::string failedLinesTmpFileNameRv = failedLines.second.substr(0, dotPos) + std::to_string(i) + failedLines.second.substr(dotPos);
        
        outFileFailedLineRV = std::make_shared<std::ofstream>(failedLinesTmpFileNameRv);
        increase_buffer_for_stream(outFileFailedLineRV);
        if (!outFileFailedLineRV->is_open()) 
        {
            std::cerr << "Error opening file: " << failedLines.second << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    std::unique_lock<std::mutex> fileLock(*threadFileOpenerMutex);
    --(*threadToInitializePtr);
    //add the new list of temporary streams to map
    tmpStreamMap[boost::this_thread::get_id()] = tmpFileStreams;

    //add the temporary failedLine to map
    failedStreamMap[boost::this_thread::get_id()] = std::make_pair(outFileFailedLineFW, outFileFailedLineRV);

    //add thread-specific statistics
    statisticsMap[boost::this_thread::get_id()] = std::make_shared<DemultiplexingStats>();

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

void DemultiplexedResult::initialize_statistiscs( const MultipleBarcodePatternVectorPtr& barcodePatternList)
{
    for(const auto& threadStatistics : statisticsMap) 
    {
        std::shared_ptr<DemultiplexingStats> stats = threadStatistics.second;
        stats->initializeStats(barcodePatternList);
    }
}

void DemultiplexedResult::initialize_thread_streams(boost::asio::thread_pool& pool, const int threadNum)
{
    //create vector that stores all custom sized (1MB) buffers for thread streams
    bufferVector.resize(4*threadNum);//buffer for threads * (dnaLine, barcodeLine, failedFilesLineFW + RV)

    //add the initialiozation fucntion thread times, where every thread waits until threadsStarted is equal
    //the number of threads to initialize every thread exactly once
    for(int i = 0; i < threadNum; ++i)
    {
        boost::asio::post(pool, std::bind(&DemultiplexedResult::initialize_tmp_file, this, i));
    }

    //only continue once all htreads have finished 
    // (we can NOT join, to not change the thread_ids, and waiting within thread only might let threadToInitialize go out of scope)
    std::unique_lock<std::mutex> lock(*threadWaitingMutex);
    cvPtr->wait(lock, [&] { return (*threadToInitializePtr) == 0; });
    lock.unlock();
}

/// write ONE failed line into a txt file
void DemultiplexedResult::write_failed_line(std::pair<std::shared_ptr<std::ofstream>, std::shared_ptr<std::ofstream>>& failedFileStream, const std::pair<fastqLine, fastqLine>& failedLine)
{
    if(failedFileStream.second == nullptr)
    {
        //for single-read write only first entry into first stream (second is a nullptr)
        *failedFileStream.first << failedLine.first.line << "\n";
    }
    else
    {
        //for paired-end sequencing write both reads
        *failedFileStream.first << failedLine.first.line << "\n";
        *failedFileStream.second << failedLine.second.line << "\n";
    }
}

/// write failed lines into a txt file
void DemultiplexedResult::write_failed_lines(std::pair<std::shared_ptr<std::ofstream>, std::shared_ptr<std::ofstream>>& failedFileStream, const std::vector<std::pair<fastqLine, fastqLine>>& failedLineBatch)
{
    std::ostringstream failedLinesFwBuffer; //temporary buffer for the batch of fw failed lines
    std::ostringstream failedLinesRvBuffer; //temporary buffer for the batch of rv failed lines

    for(const std::pair<fastqLine, fastqLine>& failedLine : failedLineBatch)
    {
        if(failedFileStream.second == nullptr)
        {
            //for single-read write only first entry into first stream (second is a nullptr)
            failedLinesFwBuffer << failedLine.first.line << "\n";
        }
        else
        {
            //for paired-end sequencing write both reads
            failedLinesFwBuffer << failedLine.first.line << "\n";
            failedLinesRvBuffer << failedLine.second.line << "\n";
        }
    }
    *failedFileStream.first << failedLinesFwBuffer.str();
    if(failedFileStream.second != nullptr)
    {
        *failedFileStream.second << failedLinesRvBuffer.str();
    }
}

/// write mapped barcodes to a tab separated file
/*
void DemultiplexedResult::write_demultiplexed_barcodes(const input& input, std::vector<std::vector<std::string>> barcodes, const std::string& patternName)
{
    std::string output = input.outPath;
    std::ofstream outputFile;

    std::string demultiplexedBarcodesFileName = patternName + ".tsv";
    if(input.prefix != "")
    {
        demultiplexedBarcodesFileName = input.prefix + "_" + demultiplexedBarcodesFileName;
    }
    std::string demultiplexedBarcodesOutput = output + "/" + demultiplexedBarcodesFileName;

    //write the barcodes we mapped
    outputFile.open (demultiplexedBarcodesOutput, std::ofstream::app);
    for(size_t i = 0; i < barcodes.size(); ++i)
    {
        for(size_t j = 0; j < barcodes.at(i).size(); ++j)
        {
            outputFile << barcodes.at(i).at(j);
            if(j!=barcodes.at(i).size()-1){outputFile << "\t";}
        }
        outputFile << "\n";
    }
    outputFile.close();
}
*/