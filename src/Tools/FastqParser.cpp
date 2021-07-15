#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <regex>
#include <thread>
#include "Barcode.hpp"

/** @param:
 *  sequence pattern: could be a sequence of constant and variable nucleotides: e.g.: [xxxxx][AGCGTACGCGAGT][xxxxx][AAGCGtAGCTTC][xxxxx] 
 *  mismatches per sequence in bracket: e.g.: 1,2,1,2,1
 * 
 * 
 **/

/*
HOW TO ALGORITHM:
    iterate over patterns and match it to substring plus/minus mismatches on both sides
    allow mismatches at beginning and end (plus/minus those mismatches, bcs imagine a barcode match with a mismatch in the end,
    we then donnt know if its  really a msimatch, or a deletion of the barcode and already part of the next barcode...)
    each matched barcodes is described by the first and last match of the sequence (therefore can be shorter, than real sequence, but not longer)
    UMI or WildcardBarcodes are matched according to the two last matches in the neighboring sequences

*/

//TODO:
/*
make a basetype of barcode:
 - derived classes of constant codes with ONE string and variables one with X strings
 - iterate through barcodes and match them in increasing order

 - if everything matches within the order of mismatches reporta vector of matches

 INPUT: pattern plus mismatches and barcode sequences
 => vector of barcode pointer which inherits to constant and variable (difference std::stirng anf vector of strings)

 OUTPUT: three files
    1.) a files with the real data strings that were amtched to barcodes
    2.) a file with the actual REAL barcode without mismatches that was mapped on the data
    3.) a file with the statistics of all barcodes, the number of mismatches per barcode

PLAN:
- design classes plus a function to calc dist for one or x sequences and return best match
- parse pattern, mismatches, and barcodeList and store in new Baseclasspointer vector

*/

#include "seqtk/kseq.h"

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

KSEQ_INIT(gzFile, gzread)

//old data types
typedef std::vector< std::shared_ptr<std::string> > BarcodeMapping;
typedef std::vector<BarcodeMapping> BarcodeMappingVector;
using namespace boost::program_options;

bool parse_arguments(char** argv, int argc, input& input)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input,i", value<std::string>(&(input.inFile))->required(), "directory of files or single file in fastq(.gz) format")
            ("output,o", value<std::string>(&(input.outFile))->required(), "output file with all split barcodes")
            
            ("sequencePattern,p", value<std::string>(&(input.patternLine))->required(), "pattern for the sequence to match, \
            every substring that should be matched is enclosed with square brackets. N is a wild card match, a substring \
            should not be a combination of wild card chars and constant chars : [AGCTATCACGTAGC][NNNNNN][AGAGCATGCCTTCAG][NNNNNN]")
            ("barcodeList,b", value<std::string>(&(input.barcodeFile)), "file with a list of all allowed well barcodes (comma seperated barcodes across several rows)\
            the row refers to the correponding bracket enclosed sequence substring. E.g. for two bracket enclosed substrings in out sequence a possible list could be:\
            AGCTTCGAG,ACGTTCAGG\nACGTCTAGACT,ATCGGCATACG,ATCGCGATC,ATCGCGCATAC")
            ("mismatches,m", value<std::string>(&(input.mismatchLine))->default_value("1"), "list of mismatches allowed for each bracket enclosed sequence substring. \
            This should be a comma seperated list of numbers for each substring of the sequence enclosed in squared brackets. E.g.: 2,1,2,1,2")
           
            ("threat,t", value<int>(&(input.threads))->default_value(5), "number of threads")

            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";
            std::cout << "EXAMPLE CALL:\n ./bin/parser -i ./inFile -o ./outFile -p [AGTCAGTC][NNNN] -m 2,1 -t 5\n";
            return false;
        }

        notify(vm);
    }
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    return true;

}

void write_file(std::string output, BarcodeMappingVector barcodes, BarcodeMappingVector realBarcodes, const std::vector<std::pair<std::string, char> > patterns)
{
    //write actually found barcodes
    std::ofstream outputFile;
    outputFile.open (output);
    //write header line
    for(int i =0; i < patterns.size(); ++i)
    {
        outputFile << patterns.at(i).first << "\t";
    }
    outputFile << "\n";

    //write all information lines
    for(int i = 0; i < barcodes.size(); ++i)
    {
        for(int j = 0; j < patterns.size(); ++j)
        {
            outputFile << *(barcodes.at(i).at(j)) << "\t";
        }
        outputFile << "\n";
    }
    outputFile.close();

    //write the barcodes we mapped
    std::size_t found = output.find_last_of("/");
    if(found == std::string::npos)
    {
        output = "CORRECTED" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "CORRECTED" + output.substr(found+1);
    }
    outputFile.open (output);
    //write header line
    for(int i =0; i < patterns.size(); ++i)
    {
        outputFile << patterns.at(i).first << "\t";
    }
    outputFile << "\n";

    //write all information lines
    for(int i = 0; i < realBarcodes.size(); ++i)
    {
        for(int j = 0; j < patterns.size(); ++j)
        {
            outputFile << *(realBarcodes.at(i).at(j)) << "\t";
        }
        outputFile << "\n";
    }
    outputFile.close();
}

void writeStats(std::string output, const fastqStats& stats)
{
    std::ofstream outputFile;
    std::size_t found = output.find_last_of("/");
    if(found == std::string::npos)
    {
        output = "STATS" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "STATS" + output.substr(found+1);
    }
    outputFile.open (output);

    for(std::pair<std::string, std::vector<int> > mismatchDictEntry : stats.mapping_dict)
    {
        outputFile << mismatchDictEntry.first << "\t";
        for(int mismatchCountIdx = 0; mismatchCountIdx < mismatchDictEntry.second.size(); ++mismatchCountIdx)
        {
            outputFile << mismatchDictEntry.second.at(mismatchCountIdx) << "\t";
        }
        outputFile << "\n";
    }

    outputFile.close();
}

//apply all BarcodePatterns to a single line to generate a vector strings (the Barcodemapping)
bool split_line_into_barcode_mappings(const std::string& seq, input* input, BarcodeMapping& barcodeMap, BarcodeMapping& realBarcodeMap,
                                      fastqStats& stats, std::map<std::string, std::shared_ptr<std::string> >& unique_seq,
                                      BarcodePatternVectorPtr barcodePatterns)
{
    //temporary vecotr of all matches
    std::vector<std::string> barcodeSequences;

    //iterate over BarcodeMappingVector
    int offset = 0;
    int score_sum = 0;

    int old_offset = offset;
    bool wildCardToFill = false;
    int wildCardLength = 0;
    for(BarcodePatternVector::iterator patternItr = barcodePatterns->begin(); 
        patternItr < barcodePatterns->end(); 
        ++patternItr)
    {
        //if we have a wildcard skip this matching, we match again the next sequence
        if((*patternItr)->is_wildcard())
        {
            old_offset = offset;
            wildCardLength = (*patternItr)->get_patterns().at(0).length();
            offset += wildCardLength;
            wildCardToFill = true;
            continue;
        }
        //for every barcodeMapping element find a match
        std::string realBarcode = ""; //the actual real barcode that we find (mismatch corrected)
        int start=0, end=0, score = 0;
        if(!(*patternItr)->match_pattern(seq, offset, start, end, score, realBarcode, stats))
        {
            std::cout << "ERROR IN: " << seq << " : " << (*patternItr)->get_patterns().at(0) << " score: " << score << "from "<< offset << "\n";
            ++stats.noMatches;
            return false;
        }
        std::cout << "pattern matched" << realBarcode << "\n";
        std::string mappedBarcode = seq.substr(offset + start, end-start);
        offset += end;
        score_sum += score;

        assert(realBarcode != "");
        //add barcode data to statistics dictionary
        int dictvectorIndex = (score <= (*patternItr)->mismatches ? score : (score + 1) );
        ++stats.mapping_dict[realBarcode].at(dictvectorIndex);
        
        //squeeze in the last wildcard match if there was one 
        //(we needed both matches of neighboring barcodes to define wildcard boundaries)
        if(wildCardToFill == true)
        {
            wildCardToFill = false;
            std::string oldWildcardMappedBarcode = seq.substr(old_offset, (offset+start-end) - old_offset);
            barcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
            realBarcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
        }
        //add this match to the BarcodeMapping
        barcodeMap.emplace_back(std::make_shared<std::string>(mappedBarcode));
        realBarcodeMap.emplace_back(std::make_shared<std::string>(realBarcode));

    }
    //if the last barcode was a WIldcard that still has to be added
    if(wildCardToFill == true)
    {
        wildCardToFill = false; // unnecessary, still left to explicitely set to false
        std::string oldWildcardMappedBarcode = seq.substr(old_offset, seq.length() - old_offset);
        barcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
        realBarcodeMap.emplace_back(std::make_shared<std::string>(oldWildcardMappedBarcode));
    }

    if(score_sum == 0)
    {
        ++stats.perfectMatches;
    }
    else
    {
        ++stats.moderateMatches;
    }
    return true;
}

void map_pattern_to_fastq_lines(const std::vector<std::string>& fastqLines, input* input, BarcodeMappingVector& barcodes, 
                                BarcodeMappingVector& realBarcodes, fastqStats& stats, BarcodePatternVectorPtr barcodePatterns)
{
    std::map<std::string, std::shared_ptr<std::string> > unique_seq;
    for(const std::string& line : fastqLines)
    {
        BarcodeMapping barcode;
        BarcodeMapping realBarcode;
        if(split_line_into_barcode_mappings(line, input, barcode, realBarcode, stats, unique_seq, barcodePatterns))
        {
            barcodes.push_back(barcode);
            realBarcodes.push_back(realBarcode);
            barcode.clear();
        }
    }
}

void iterate_over_fastq(input& input, BarcodeMappingVector& barcodes, BarcodeMappingVector& realBarcodes, BarcodePatternVectorPtr barcodePatterns, fastqStats& fastqStatsFinal)
{
    //read all fastq lines into str vector
    gzFile fp;
    kseq_t *ks;
    fp = gzopen(input.inFile.c_str(),"r");
    if(NULL == fp){
        fprintf(stderr,"Fail to open file: %s\n", input.inFile.c_str());
    }
    ks = kseq_init(fp);
    std::vector<std::string> fastqLines;
    while( kseq_read(ks) >= 0 ){
        fastqLines.push_back(std::string(ks->seq.s));
    }

    //split input lines into thread buckets
    std::vector<std::vector<std::string> > fastqLinesVector; //vector holding all the buckets of fastq lines to be analyzed by each thread
                                                            // each vector element is one bucket for one thread
    int element_number = fastqLines.size() / input.threads;
    std::vector<std::string>::iterator begin = fastqLines.begin();
    for(int i = 0; i < input.threads; ++i)
    {
        std::vector<std::string>::iterator end = ( i==input.threads-1 ? fastqLines.end() : begin + element_number);

        std::vector<std::string> tmpFastqLines(begin, end);
        fastqLinesVector.push_back(tmpFastqLines);
        begin = end;
    }

    //for every bucket call a thread
    //tmp variables to store thread results
    std::vector<BarcodeMappingVector> barcodesThreadList(input.threads); // vector of mappedSequences to be filled
    std::vector<BarcodeMappingVector> realBarcodesThreadList(input.threads); // vector of mappedSequences to be filled with corrected matches
    std::vector<fastqStats> statsThreadList(input.threads, fastqStatsFinal); // vector of stats to be filled
    std::vector<std::thread> workers;
    for (int i = 0; i < input.threads; ++i) {
        workers.push_back(std::thread(map_pattern_to_fastq_lines, fastqLinesVector.at(i), 
                        &input, std::ref(barcodesThreadList.at(i)), std::ref(realBarcodesThreadList.at(i)), 
                        std::ref(statsThreadList.at(i)), barcodePatterns));
    }
    for (std::thread &t: workers) 
    {
        if (t.joinable()) {
            t.join();
        }
    }
    //combine thread data
    for (int i = 0; i < input.threads; ++i) 
    {
        barcodes.insert(barcodes.end(), barcodesThreadList.at(i).begin(), barcodesThreadList.at(i).end());
        realBarcodes.insert(realBarcodes.end(), realBarcodesThreadList.at(i).begin(), realBarcodesThreadList.at(i).end());
        fastqStatsFinal.perfectMatches += statsThreadList.at(i).perfectMatches;
        fastqStatsFinal.noMatches += statsThreadList.at(i).noMatches;
        fastqStatsFinal.moderateMatches += statsThreadList.at(i).moderateMatches;
        fastqStatsFinal.multiBarcodeMatch += statsThreadList.at(i).multiBarcodeMatch;
        //iterate for every thread over the dictionary of mismatches per barcode and fill the final stats data with it
        for(std::pair<std::string, std::vector<int> > mismatchDictEntry : statsThreadList.at(i).mapping_dict)
        {
            for(int mismatchCountIdx = 0; mismatchCountIdx < mismatchDictEntry.second.size(); ++mismatchCountIdx)
            {
                fastqStatsFinal.mapping_dict[mismatchDictEntry.first].at(mismatchCountIdx) += statsThreadList.at(i).mapping_dict[mismatchDictEntry.first].at(mismatchCountIdx);
            }
        }
    }
    std::cout << "MATCHED: " << fastqStatsFinal.perfectMatches << " | MODERATE MATCH: " << fastqStatsFinal.moderateMatches
              << " | MISMATCHED: " << fastqStatsFinal.noMatches << " | Multiplebarcode: " << fastqStatsFinal.multiBarcodeMatch << "\n";
    kseq_destroy(ks);
    gzclose(fp);
}

void parseBarcodeData(const input& input, std::vector<std::pair<std::string, char> >& patterns, std::vector<int>& mismatches, std::vector<std::vector<std::string> >& varyingBarcodes)
{
    try{
        // parse the pattern, mismatches, and barcode file (perform quality check as well)
        int numberOfNonConstantBarcodes = 0;
        std::string pattern = input.patternLine;
        std::string delimiter = "]";
        const char delimiter2 = '[';
        size_t pos = 0;
        std::string seq;
        //PARSE PATTERN
        std::pair<std::string, bool> seqPair;
        while ((pos = pattern.find(delimiter)) != std::string::npos) {
            seq = pattern.substr(0, pos);
            pattern.erase(0, pos + 1);
            if(seq.at(0) != delimiter2)
            {
                std::cerr << "PARAMETER ERROR: Wrong barcode patter parameter, it looks like a \'[\' is missing\n";
                exit(1);
            }
            seq.erase(0, 1);
            //just check if we have a barcode pattern of only'N', bcs the nunber of those patterns must match the number of lines in barcodefile
            bool nonConstantSeq = true;
            char patternType = 'c';
            for (char const &c: seq) {
                if(!(c=='A' | c=='T' | c=='G' |c=='C' |
                    c=='a' | c=='t' | c=='g' | c=='c' |
                    c=='N' | c=='X'))
                {
                    std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C,N)\n";
                    if(c==' ' | c=='\t' | c=='\n')
                    {
                        std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
                    }
                    exit(1);
                }
                if(c!='N'){nonConstantSeq=false;}

                //determine pattern type and throw error
                if( (patternType!='c') & (c!='N') & (c!='X') )
                {
                    std::cerr << "PARAMETER ERROR: a barcode sequence has bases as well as wildcard 'X' or variable 'N' bases, the combination is not allowed!!!\n";
                    exit(1);
                }
                if(c=='N'){patternType='v';}
                else if(c=='X'){patternType='w';}
            }
            if(nonConstantSeq){++numberOfNonConstantBarcodes;}
            
            patterns.push_back(std::make_pair(seq, patternType));
        }
        //PARSE mismatches
        pattern = input.mismatchLine;
        delimiter = ",";
        pos = 0;
        while ((pos = pattern.find(delimiter)) != std::string::npos) {
            seq = pattern.substr(0, pos);
            pattern.erase(0, pos + 1);
            mismatches.push_back(stoi(seq));
        }
        mismatches.push_back(stoi(pattern));
        if(patterns.size() != mismatches.size())
        {
            std::cerr << "PARAMETER ERROR: Number of barcode patterns and mismatches is not equal\n";
            exit(1);
        }
        //PARSE barcode file
        std::ifstream barcodeFile(input.barcodeFile);
        for(std::string line; std::getline(barcodeFile, line);)
        {
            delimiter = ",";
            pos = 0;
            std::vector<std::string> seqVector;
            while ((pos = line.find(delimiter)) != std::string::npos) {
                seq = line.substr(0, pos);
                line.erase(0, pos + 1);
                for (char const &c: seq) {
                    if(!(c=='A' | c=='T' | c=='G' |c=='C' |
                         c=='a' | c=='t' | c=='g' | c=='c'))
                         {
                            std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C)\n";
                            if(c==' ' | c=='\t' | c=='\n')
                            {
                                std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
                            }
                            exit(1);
                         }
                }
                seqVector.push_back(seq);
            }
            seq = line;
            for (char const &c: seq) {
                if(!(c=='A' | c=='T' | c=='G' |c=='C' |
                        c=='a' | c=='t' | c=='g' | c=='c'))
                        {
                        std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C)\n";
                        if(c==' ' | c=='\t' | c=='\n')
                        {
                            std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
                        }
                        exit(1);
                        }
            }
            seqVector.push_back(seq);
            varyingBarcodes.push_back(seqVector);
            seqVector.clear();
        }
        if(numberOfNonConstantBarcodes != varyingBarcodes.size())
        {
            std::cerr << "PARAMETER ERROR: Number of barcode patterns for non-constant sequences [N*] and lines in barcode file are not equal\n";
            exit(1);
        }
        
    }
    catch(std::exception& e)
    {
        std::cerr << "PARAMETER ERROR: Check the input parameters again (barcode pattern list, mismatch list, file with barcode lists)\n";
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

BarcodePatternVectorPtr generate_barcode_patterns(input input, std::vector<std::pair<std::string, char> >& patterns)
{
    std::vector<int> mismatches; // vector of all string patterns
    std::vector<std::vector<std::string> > varyingBarcodes; // a vector storing for all non-constant barcode patterns in the order of occurence
                                                            // in the barcode pattern string the possible barcode sequences
    //fill the upper three vectors and handle as many errors as possible
    parseBarcodeData(input, patterns, mismatches, varyingBarcodes);

/*
//DEBUG OUTPUT TO CHECK PARSED PARAMETERS TO BUILD BARCODE PATTERNS
    int x = 0;
    for(int i=0; i< patterns.size(); ++i)
    {
        std::cout << patterns.at(i).first << "_"<< patterns.at(i).second << "\n";
        std::cout << mismatches.at(i) << "\n";
        if(patterns.at(i).second)
        {
            for(int j =0; j<varyingBarcodes.at(x).size(); ++j)
            {
                std::cout << varyingBarcodes.at(x).at(j);
            }
            std::cout << "\n";
            ++x;
        }
    }
*/
    //iterate over patterns and fill BarcodePatternVector instance
    BarcodePatternVector barcodeVector;
    int variableBarcodeIdx = 0; // index for variable barcode sequences is shorter than the vector of patterns in total
    for(int i=0; i < patterns.size(); ++i)
    {
        if(patterns.at(i).second=='v')
        {
            VariableBarcode barcode(varyingBarcodes.at(variableBarcodeIdx), mismatches.at(i));
            std::shared_ptr<VariableBarcode> barcodePtr(std::make_shared<VariableBarcode>(barcode));
            barcodeVector.emplace_back(barcodePtr);
            ++variableBarcodeIdx;
        }
        else if(patterns.at(i).second=='c')
        {
            ConstantBarcode barcode(patterns.at(i).first, mismatches.at(i));
            std::shared_ptr<ConstantBarcode> barcodePtr(std::make_shared<ConstantBarcode>(barcode));
            barcodeVector.emplace_back(barcodePtr);
        }
        else if(patterns.at(i).second=='w')
        {
            WildcardBarcode barcode(patterns.at(i).first, mismatches.at(i));
            std::shared_ptr<WildcardBarcode> barcodePtr(std::make_shared<WildcardBarcode>(barcode));
            barcodeVector.emplace_back(barcodePtr);
        }
    }
    BarcodePatternVectorPtr barcodePatternVector = std::make_shared<BarcodePatternVector>(barcodeVector);

    return barcodePatternVector;
}

//TODO: so far only implemented for a single fastq file,
//implement a function to iterate over all fastq files in a directory and combine data
//e.g. if we have several files for different batches etc. ...
void split_barcodes(input input, BarcodeMappingVector& barcodes, BarcodeMappingVector& realBarcodes, BarcodePatternVectorPtr barcodePatterns, fastqStats& fastqStatsFinal)
{
    //if directory iterate over all fastqs

        //combine data to one

    //if file
    iterate_over_fastq(input, barcodes, realBarcodes, barcodePatterns, fastqStatsFinal);

}

//initialize the dictionary of mismatches in each specific barcode
// the vector in the dict has length "mismatches + 2" for each barcode
// one entry for zero mismatches, eveery number from, 1 to mismatches and 
//one for more mismatches than the max allowed number of mismathces
void initializeStats(fastqStats& stats, const BarcodePatternVectorPtr barcodePatterns)
{
    for(BarcodePatternVector::iterator patternItr = barcodePatterns->begin(); 
        patternItr < barcodePatterns->end(); 
        ++patternItr)
    {
        int mismatches = (*patternItr)->mismatches;
        const std::vector<std::string> patterns = (*patternItr)->get_patterns();
        for(const std::string& pattern : patterns)
        {
            std::vector<int> mismatchVector(mismatches + 2, 0);
            stats.mapping_dict.insert(std::make_pair(pattern, mismatchVector));
        }
    }
}

int main(int argc, char** argv)
{
    input input;
    std::vector<std::pair<std::string, char> > patterns; // vector of all string patterns, 
                                                        //second entry is c=constant, v=varying, w=wildcard
    BarcodeMappingVector mappedBarcodes;
    BarcodeMappingVector realBarcodes;
    fastqStats fastqStats;
    if(parse_arguments(argv, argc, input))
    {
        BarcodePatternVectorPtr barcodePatterns = generate_barcode_patterns(input, patterns);
        initializeStats(fastqStats, barcodePatterns);
        split_barcodes(input, mappedBarcodes, realBarcodes, barcodePatterns, fastqStats);
        write_file(input.outFile, mappedBarcodes, realBarcodes, patterns);
        writeStats(input.outFile, fastqStats);
    }
 
    return EXIT_SUCCESS;
}