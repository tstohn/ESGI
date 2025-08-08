#pragma once

#include "DemultiplexedResult.hpp"
#include <limits>

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

        /*void demultiplex_wrapper(const std::pair<fastqLine, fastqLine>& line,
                                const input& input,
                                std::atomic<int>& lineCount,
                                const unsigned long long& totalReadCount,
                                std::atomic<long long int>& elementsInQueue);*/
        void demultiplex_wrapper_batch(const std::vector<std::pair<fastqLine, fastqLine>>& line_vector,
                                                const input& input,
                                                std::atomic<unsigned long long>& lineCount,
                                                const unsigned long long& totalReadCount,
                                                std::atomic<long long int>& elementsInQueue);
        void run_mapping(const input& input);

        //map to store the temporary output files (e.g., for RAM efficient laptop usage)
        //theadID maps to a vector of ordered fileStreams for every pattern in the order of
        //the inout file. ASSURE THE RIGHT ORDER!! We explicitely do not use the patternNames
        //to avoid another string-hash as we can also simply keep the order of patterns
        DemultiplexedResultPtr fileWriter;

        std::unordered_map<boost::thread::id, MultipleBarcodePatternVectorPtr, thread_id_hash> thread_pattern;
        
        void create_pattern_copy_in_thread(std::shared_ptr<std::mutex> threadFillMutex,
                                           std::shared_ptr<std::mutex> threadWaitingMutex,
                                           std::shared_ptr<std::atomic<unsigned int> > threadToInitializePtr,
                                           std::shared_ptr<std::condition_variable> cvPtr)
        {
            MultipleBarcodePatternVectorPtr origional = this->get_barcode_pattern();
            MultipleBarcodePatternVectorPtr copy = std::make_shared<std::vector<BarcodePatternPtr>>();

            for (const auto& pattern : *origional)
            {
                BarcodePatternPtr patternCopyPtr = std::make_shared<BarcodePattern>(*pattern);
                copy->push_back(pattern ? patternCopyPtr : nullptr);
            }

            std::unique_lock<std::mutex> threadFillLock(*threadFillMutex);
            --(*threadToInitializePtr);
            //add the pattern copy to thead->copy map
            thread_pattern.emplace(boost::this_thread::get_id(), copy);

            //decrease number of threads that need initialization, when all are initialized we can continue program in main function
            if (*threadToInitializePtr == 0) 
            {
                cvPtr->notify_all();  // Notify when all tasks are done
                threadFillLock.unlock();
            }
            else
            {
                //lock threads if we did not finish them - like this we can be 100% sure that evey thread gets called EXACTLY ONCE
                // for initialization
                threadFillLock.unlock();
                //wait until all threads are here
                std::unique_lock<std::mutex> lock(*threadWaitingMutex);
                cvPtr->wait(lock, [&] { return (*threadToInitializePtr) == 0; });
                lock.unlock();
            }
        }

        void initialize_thread_patterns(boost::asio::thread_pool& pool, const int threadNum)
        {
            std::shared_ptr<std::mutex> threadFillMutex = std::make_shared<std::mutex>();
            std::shared_ptr<std::mutex> threadWaitingMutex = std::make_shared<std::mutex>();
            std::shared_ptr<std::atomic<unsigned int>> threadToInitializePtr = std::make_shared<std::atomic<unsigned int>>(threadNum);
            std::shared_ptr<std::condition_variable> cvPtr = std::make_shared<std::condition_variable>();
  
            for(int i = 0; i < threadNum; ++i)
            {
                boost::asio::post(pool, std::bind(&Demultiplexer::create_pattern_copy_in_thread, this, threadFillMutex, threadWaitingMutex, 
                threadToInitializePtr, cvPtr));
            }

            //only continue once all htreads have finished 
            // (we can NOT join, to not change the thread_ids, and waiting within thread only might let threadToInitialize go out of scope)
            std::unique_lock<std::mutex> lock(*threadWaitingMutex);
            cvPtr->wait(lock, [&] { return (*threadToInitializePtr) == 0; });
            lock.unlock();
        }

    public:
        void run(const input& input);

};