#include "BarcodeMapping.hpp"

/** @brief class overriting a couple of functions of Mapping class 
 * to store statistics, failes lines, etc
 * also this class allows to read only a subset of reads into RAM
 * and by that keeping the processing queue only filled up to a certain
 * level, once a task is finished the queue is filled with the next read
**/
template<typename MappingPolicy, typename FilePolicy>
class DemultiplexedLinesWriter : private Mapping<MappingPolicy, FilePolicy>
{
    private:

        void demultiplex_wrapper(const std::string& line,
                                const input& input,
                                std::atomic<unsigned long long>& lineCount,
                                const unsigned long long& totalReadCount,
                                std::atomic<long long int>& elementsInQueue);
        void initialize_output_files(const input& input,const std::vector<std::pair<std::string, char> >& patterns);
        void run_mapping(const input& input);


    public:
        void run(const input& input);

};