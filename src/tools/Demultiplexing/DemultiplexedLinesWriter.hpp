#include "mapping.hpp"

//class writing the demultiplexing data to files
//handles smaller chunks of data that are written to a file immediately, 
//and stores additional information as: lines the failed demultipelxing,
//the real sequence to which a barcode was mapped, etc.

template<typename MappingPolicy, typename FilePolicy>
class DemultiplexedLinesWriter : private Mapping<MappingPolicy, FilePolicy>
{
    private:

        std::shared_ptr<fastqStats> fastqStatsPtr;

        void initialize_mapping(const input& input,const std::vector<std::pair<std::string, char> >& patterns);
        void run_mapping(const input& input);


    public:
        void run(const input& input);

};