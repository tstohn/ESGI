#pragma once

#include "OutputFileWriter.hpp"

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