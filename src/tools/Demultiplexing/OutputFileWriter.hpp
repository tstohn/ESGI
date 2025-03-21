#pragma once

#include "BarcodeMapping.hpp"

#include <condition_variable>
#include <cstdio>  // For std::remove()

//threadID hash to map every thread to its own temporary files to write
struct thread_id_hash {
    size_t operator()(const boost::thread::id& id) const {
        return boost::hash<boost::thread::id>()(id);
    }
  };
  
  //a pattern can have a barcode & maybe a dna region
  // therefore every thread must be able to potentially write to both files 
  // this struct saves both and contains a nullptr in case of absence
  struct TmpPatternStream
  {
      std::shared_ptr<std::ofstream> barcodeStream = nullptr;
      std::shared_ptr<std::ofstream> dnaStream = nullptr;
      unsigned long lineNumber = 0;
  };
  
  struct FinalPatternFiles
  {
      std::string barcodeFile = "";
      std::string dnaFile = "";
  };
  /** @brief class to handle thread-specific output streams
  * streams are initialized from pattern data
  * depending on weather pattern contains only barcodes or also dna
  * the threadId can provide a barcode & dna stream to write to
  * finally this class can concatenate final files
  */
  // ONE FAiledFile per thread -> in the end combined into one
  // for every thread a map of pattern (which contains DNA) to its DNA/BARCODE file
  class OutputFileWriter
  {
      public:
  
          //the final files are created upon initilization
          //HOWEVER, tmp files per thread need to be created when the thread-pool is set up
          OutputFileWriter(const input& input, const MultipleBarcodePatternVectorPtr& barcodePatternList)
          {
              //initialze the names/ headers of final output files for each pattern
              initialize(input, barcodePatternList);
              //thread-dependent tmp files are initialized later
              //but initialyze the mutax for the htread initialization, which has to be shared for copy-construction of the outputFilewriter
              threadFileOpenerMutex = std::make_unique<std::mutex>();
              threadWaitingMutex = std::make_unique<std::mutex>();
  
              threadToInitializePtr = std::make_unique<std::atomic<unsigned int>>(input.threads);
              cvPtr = std::make_unique<std::condition_variable>();
  
          }
          void initialize_thread_streams(boost::asio::thread_pool& pool, const int threadNum);
  
          //return the tmp-stream for barcodes of the pattern-name for a threadID
          std::shared_ptr<std::ofstream> get_barcodeStream_for_threadID_at(const boost::thread::id& threadID, const std::string& patternName)
          {
              return(tmpStreamMap.at(threadID).at(patternName).barcodeStream);
          }
          //return the tmp-stream for DNA of the pattern-name for a threadID
          std::shared_ptr<std::ofstream> get_dnaStream_for_threadID_at(const boost::thread::id& threadID, const std::string& patternName)
          {
              return(tmpStreamMap.at(threadID).at(patternName).dnaStream);
          }
        //return the tmp-stream for DNA of the pattern-name for a threadID
        TmpPatternStream& get_streams_for_threadID(const boost::thread::id& threadID, const std::string& patternName)
        {
            return(tmpStreamMap.at(threadID).at(patternName));
        }

          //return the tmp-stream for DNA of the pattern-name for a threadID
          std::pair<std::shared_ptr<std::ofstream>, std::shared_ptr<std::ofstream>> get_failedStream_for_threadID_at(const boost::thread::id& threadID)
          {
              return(failedStreamMap.at(threadID));
          }

          //return the demultiplexed reads for a pattern
          const BarcodeMappingVector get_demultiplexed_ab_reads(const std::string& patternName)
          {
              return(demultiplexedReads.at(patternName)->get_all_reads());
          }
          void add_demultiplexed_line(const std::string& patternName, std::vector<std::string>& demultiplexedLineString)
          {
            demultiplexedReads.at(patternName)->addVector(demultiplexedLineString);
          }
  
          void concatenateFiles(const std::vector<std::string>& tmpFileList, const std::string& outputFile);
          //writing of final files
          void close_and_concatenate_fileStreams(const input& input);
  
          //write a failed line that is encountered to the tmp-threadFileStreams
          void write_failed_line(std::pair<std::shared_ptr<std::ofstream>, std::shared_ptr<std::ofstream>>& failedFileStream, const std::pair<fastqLine, fastqLine>& failedLine);
          
          //write the fastq and barcode data into temporary stream
          void write_dna_line(TmpPatternStream& dnaLineStream, const DemultiplexedLine& demultiplexedLine, const boost::thread::id& threadID);

          //write final files: from memory or by concatenating & deleting tmp-files
          void write_demultiplexed_barcodes(const input& input, BarcodeMappingVector barcodes, const std::string& patternName);

          // 1.) Write DemutltiplexedReads as barcode-files for every pattern
          // 2.) write the 2 mismatch files: a) mismatches per barcode b) mismatches for the different patterns
          // 3.) close thw failedLines streams per thread, write into ONE output, and delte tmp files (file per thread)
          // 4.) for every thread, fo every pattern (with DNA) close the tmp-files, write to a final
          // DNa(fastq)&barcode(TSV) file, and delete tmp files (file per thread per dna-pattern)
          void write_output(const input& input);
  
      private:
  
          void initialize_additional_output(const input& input);
          void initialize_output_for_pattern(std::string output, const BarcodePatternPtr pattern);
          //initializes the files for output
          void initialize(const input& input, const MultipleBarcodePatternVectorPtr& barcodePatternList);
          void initialize_tmp_file(const int i);
          
          //maps a patternName to a list of all demultipelx-reads found for this pattern
          std::unordered_map<std::string, DemultiplexedReadsPtr> demultiplexedReads;
  
          //names of final files
          std::string barcodeMismatches; //stored and written in the end
          std::string patternMismatches; // stored and written in the end
  
          //map of patternName to FinalPattern file struct 
          // the struct stores the names of the files per pattern: a demultiplexed barcode file or
          // (if DNA is contained) a fastq and a TSV file (TSV storing the barcodes with read name) 
          //the string is the same as pattern->patternName
          std::unordered_map<std::string, FinalPatternFiles> finalFiles;
          std::pair<std::string, std::string>  failedLines; //fw and rv read of failed lines, if not paired end the failed lines are stored in first entry of pair
  
          //maps threadID -> list of streams for all patterns (ordered)
          std::unordered_map<boost::thread::id, std::unordered_map<std::string, TmpPatternStream>, thread_id_hash> tmpStreamMap;
          std::unordered_map<boost::thread::id, std::pair<std::shared_ptr<std::ofstream>,std::shared_ptr<std::ofstream>>, thread_id_hash> failedStreamMap;
         
          //data structures for initialization of thread-specific temporary files
          std::unique_ptr<std::mutex> threadFileOpenerMutex;  //locking writing access to the above tmpStreamMap
          std::unique_ptr<std::mutex> threadWaitingMutex;  //locking the waiting locks to call every thread once
          //we need both variables as class variables - like this the variables exist if main thread goes on
          //while initializer threads still finish the initialization (initialization of shared variables is done at this point)
          //we only close the function and end the waiting to re-use thread for upcoming demultiplexing
          std::unique_ptr<std::atomic<unsigned int>> threadToInitializePtr; //this gets counted down in initialize_tmp_file, when zero we continue
          std::unique_ptr<std::condition_variable> cvPtr;
  
  };
  typedef std::shared_ptr<OutputFileWriter> OutputFileWriterPtr; 