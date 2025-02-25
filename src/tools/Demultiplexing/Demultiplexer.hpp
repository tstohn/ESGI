#include "BarcodeMapping.hpp"

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
    std::ofstream* barcodeStream = nullptr;
    std::ofstream* dnaStream = nullptr;
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
class OutputFileWriter
{
    public:

        OutputFileWriter(const input& input, const MultipleBarcodePatternVectorPtr& barcodePatternList)
        {
            //initialze the names/ headers of final output files for each pattern
            initialize_output_files(input, barcodePatternList);
            //thread-dependent tmp files are initialized later
            //but initialyze the mutax for the htread initialization, which has to be shared for copy-construction of the outputFilewriter
            threadFileOpenerMutex = std::make_unique<std::mutex>();
            threadWaitingMutex = std::make_unique<std::mutex>();

            threadToInitializePtr = std::make_unique<std::atomic<unsigned int>>(input.threads);
            cvPtr = std::make_unique<std::condition_variable>();

        }
        void initialize_thread_streams(boost::asio::thread_pool& pool, const int threadNum);

        //return the tmp-stream for barcodes of the pattern-name for a threadID
        std::ofstream* get_barcodeStream_for_threadID_at(const boost::thread::id& threadID, const std::string& patternName)
        {
            return(tmpStreamMap.at(threadID).at(patternName).barcodeStream);
        }
        //return the tmp-stream for DNA of the pattern-name for a threadID
        std::ofstream* get_dnaStream_for_threadID_at(const boost::thread::id& threadID, const std::string& patternName)
        {
            return(tmpStreamMap.at(threadID).at(patternName).dnaStream);
        }

        void close_and_concatenate_fileStreams()
        {
            
        }

        //write final files: from memory or by concatenating&deleting tmp-files
        void write_output()
        {
            //if tmp files were written (for failed/ DNA reads)
            close_and_concatenate_fileStreams();

            //write all barcode-only files (data is still in memory at this point)
            //write_file(input, this->get_demultiplexed_ab_reads());

            //write 2 mismatch files(pattern->MM, barcodes->MM)
            //write_stats(input, this->get_mismatch_dict());
        }

    private:

        void initialize_additional_output(std::string output);
        void initialize_output_for_pattern(std::string output, const BarcodePatternPtr pattern);
        //initializes the files for output
        void initialize_output_files(const input& input, const MultipleBarcodePatternVectorPtr& barcodePatternList);
        void initialize_tmp_file();
        
        //maps a patternName to a list of all demultipelx-reads found for this pattern
        std::unordered_map<std::string, std::vector<DemultiplexedReadsPtr>> demultiplexedReads;

        //names of final files
        std::string barcodeMismatches;
        std::string patternMismatches;
        std::string failedLines;

        //map of patternName to FinalPattern file struct (storing potential barcode and dna file)
        std::unordered_map<std::string, FinalPatternFiles> finalFiles;
        //maps threadID -> list of streams for all patterns (ordered)
        std::unordered_map<boost::thread::id, std::unordered_map<std::string, TmpPatternStream>, thread_id_hash> tmpStreamMap;
        std::unique_ptr<std::mutex> threadFileOpenerMutex;  //locking writing access to the above tmpStreamMap
        std::unique_ptr<std::mutex> threadWaitingMutex;  //locking the waiting locks to call every thread once

        //we need both variables as class variables - like this the variables exist if main thread goes on
        //while initializer threads still finish the initialization (initialization of shared variables is done at this point)
        //we only close the function and end the waiting to re-use thread for upcoming demultiplexing
        std::unique_ptr<std::atomic<unsigned int>> threadToInitializePtr; //this gets counted down in initialize_tmp_file, when zero we continue
        std::unique_ptr<std::condition_variable> cvPtr;

};
typedef std::shared_ptr<OutputFileWriter> OutputFileWriterPtr; 

/** @brief class to map several barcode Patterns simultaneously, 
 * and handles writing of results/ or storage in RAM
 * this calss is overriting a couple of functions of Mapping class 
 * to store statistics, failes lines, etc
 * also this class allows to read only a subset of reads into RAM
 * and by that keeping the processing queue only filled up to a certain
 * level, once a task is finished the queue is filled with the next read
**/
template<typename MappingPolicy, typename FilePolicy>
class Demultiplexer : private Mapping<MappingPolicy, FilePolicy>
{
    private:

        void demultiplex_wrapper(std::pair<const std::string&, const std::string&> line,
                                const input& input,
                                const unsigned long long lineCount,
                                const unsigned long long& totalReadCount,
                                std::atomic<long long int>& elementsInQueue);
        void run_mapping(const input& input);

        //map to store the temporary output files (e.g., for RAM efficient laptop usage)
        //theadID maps to a vector of ordered fileStreams for every pattern in the order of
        //the inout file. ASSURE THE RIGHT ORDER!! We explicitely do not use the patternNames
        //to avoid another string-hash as we can also simply keep the order of patterns
        OutputFileWriterPtr fileWriter;

        //dictionary mapping the barcodePattern-name to the file in which to write it in the end (final output)
        std::vector<std::string> outputFileNames;

    public:
        void run(const input& input);

};