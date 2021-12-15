#include "mapping.hpp"

//class overriting a couple of functions of mapping to store statistics,
//failes lines, etc

//class writing the demultiplexing data to files
//handles smaller chunks of data that are written to a file immediately, 
//and stores additional information (e.g. lines the failed demultipelxing)
template<typename MappingPolicy, typename FilePolicy>
class DemultiplexedLinesWriter : private Mapping<MappingPolicy, FilePolicy>
{
    private:

        void demultiplex_wrapper(const std::string& line,
                                const input& input,
                                std::atomic<int>& lineCount,
                                const int& totalReadCount,
                                std::atomic<int>& elementsInQueue);
        void initialize_output_files(const input& input,const std::vector<std::pair<std::string, char> >& patterns);
        void run_mapping(const input& input);


    public:
        void run(const input& input);

};