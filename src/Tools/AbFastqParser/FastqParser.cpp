#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <zlib.h>
#include <regex>
#include <thread>

#include "Barcode.hpp"

/** 
 * A tool to map fastq lines (stitched to one read, use e.g. fastq-join) to a certain barcode pattern,
 * this pattern must be given as an input parameter plus additional files for patterns with multiple
 * possible barcodes
 * 
 * HOW TO ALGORITHM:
    iterate over patterns and match it to substring plus/minus mismatches on both sides
    allow mismatches at beginning and end (plus/minus those mismatches, bcs imagine a barcode match with a mismatch in the end,
    we then donnt know if its  really a msimatch, or a deletion of the barcode and already part of the next barcode...)
    each matched barcodes is described by the first and last match of the sequence (therefore can be shorter, than real sequence, but not longer)
    UMI or WildcardBarcodes are matched according to the two last matches in the neighboring sequences
 * @param <input> input fastq file (gzipped), considers only full length fastq reads, FW/RV must be stitched together in advance
 * @param <output> output extension, that will be added to the output files, see return
 * @param <sequencePattern> a string with all the barcode patterns, each pattern is enclosed by suqare brackets, valid chars are AGTC ofr bases, N for a sequences
 *                          that hold a fixed number of different combinations, must be declared in the barcodeList file, and X for wildcard sequences, that can
 *                          be anything (e.g. for UMIs): e.g.: [xxxxx][AGCGTACGCGAGT][xxxxx][AAGCGtAGCTTC][xxxxx] 
 * @param <barcodeList> a file of discrete barcodes, that can be mapped to the [N...] sequences, each row is one pattern, they r in the order as in the sequence pattern and
 *                      barcodes must be comma seperated
 * @param <mismatches> comma seperated list of mismatches, one entry for each sequence pattern declared above: e.g.: 1,2,1,2,1
 *  
 * @return  writes a tab seperated file with barcodes, a file with statistical information, and if flag <-r> set, also a tab seperated file
 *          of corresponding reads with their real sequence, including mismatches.
 *          After running a line of perfect, moderate and mismatches is printed. Each count refers to a whole line (all barcodes matched perfectly,
 *          at least one matched only with mismatches, at least one did not match at all. The number of duplicate barcodes is per barcode, and thus might
 *          be greater than the total number of mismatches.)
 * 
 *          barcode file: Output file starts with 'BarcodeMapping_' and ends with the ending extension declared with the -o paramter
 *          stats file: a list where for each barcode in the fastq file, the number of found mismatches is written, the columns start with zero mismatches
 *          to the maximum number of mismtaches declared in mismatches parameter plus one, this last columns sums up all cases of more than max-num mismatches.
 *          real barcode file: like barcode file with uncorrected sequences, starts with RealBarcodeSequence
 **/

#include "seqtk/kseq.h"

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

KSEQ_INIT(gzFile, gzread)

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
            ("realSequences,r", value<bool>(&(input.storeRealSequences))->default_value(false), "set flag if next to the mapped barcodes also a table of the actual \
            reads - including mismatches - should be printed\n")
            ("fastqReadBucketSize,s", value<int>(&(input.fastqReadBucketSize))->default_value(10000000), "number of lines of the fastQ file that should be read into RAM \
            and be processed, before the next batch of fastq reads is read into RAM and processed.")

            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";

            std::cout << "###########################################\n";
            std::cout << "EXAMPLE CALL:\n ./bin/parser -i ./inFile -o ./outFile -p [AGTCAGTC][NNNN] -b ./barcodeFile.txt -m 2,1 -t 5\n";
            std::cout << "Calling the Tool, mapping each line to a pattern, that starts with a constant sequence of <AGTCAGTC> in which up to two mismatches\n\
            are allowed. After that the reads should contain a four base long sequence with a maximum of one mismatch. All sequences that are allowed are in the text file barcodeFile.txt.\n\
            This file should have in its first row, all comma seperated sequences that could map, and they should be all only four bases long, allowed chars are only A,C,G,T.\n";
            std::cout << "###########################################\n";

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

//TODO: function to calculate maximum distance between barcodes: use this distance-1 as mismatch threshold
//return value is a string of all those values for all different barcodes: b4_threshold, ab_threshold, bc1_thrshold, ...
void calculate_string_of_mismatches(std::vector<int>& mismatches, const std::string& mismatchString)
{
    
}

void write_file(const input& input, BarcodeMappingVector barcodes, BarcodeMappingVector realBarcodes, const std::vector<std::pair<std::string, char> > patterns)
{
    std::string output = input.outFile;
    std::ofstream outputFile;

    //write real sequences that map to barcodes
    std::size_t found = output.find_last_of("/");
    if(input.storeRealSequences)
    {
        std::string outputReal;
        if(found == std::string::npos)
        {
            outputReal = "RealBarcodeSequence_" + output;
        }
        else
        {
            outputReal = output.substr(0,found) + "/" + "RealBarcodeSequence_" + output.substr(found+1);
        }
        outputFile.open (outputReal, std::ofstream::app);
        for(int i = 0; i < barcodes.size(); ++i)
        {
            for(int j = 0; j < patterns.size(); ++j)
            {
                outputFile << *(barcodes.at(i).at(j));
                if(j!=patterns.size()-1){outputFile << "\t";}
            }
            outputFile << "\n";
        }
        outputFile.close();
    }
    
    //write the barcodes we mapped
    if(found == std::string::npos)
    {
        output = "BarcodeMapping_" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "BarcodeMapping_" + output.substr(found+1);
    }
    outputFile.open (output, std::ofstream::app);
    for(int i = 0; i < realBarcodes.size(); ++i)
    {
        for(int j = 0; j < patterns.size(); ++j)
        {
            outputFile << *(realBarcodes.at(i).at(j));
            if(j!=patterns.size()-1){outputFile << "\t";}
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
        output = "StatsBarcodeMappingErrors_" + output;
    }
    else
    {
        output = output.substr(0,found) + "/" + "StatsBarcodeMappingErrors_" + output.substr(found+1);
    }
    std::remove(output.c_str());
    outputFile.open (output, std::ofstream::app);

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
            ++stats.noMatches;
            return false;
        }

        std::string mappedBarcode = seq.substr(offset + start, end-start);
        offset += end;
        score_sum += score;

        assert(realBarcode != "");
        //add barcode data to statistics dictionary
        int dictvectorIndex = (score <= (*patternItr)->mismatches ? score : ( ((*patternItr)->mismatches) + 1) );
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

void writeCurrentResults(BarcodeMappingVector& barcodes, BarcodeMappingVector& realBarcodes, fastqStats& fastqStatsFinal, input& input, 
                         const std::vector<std::pair<std::string, char> >& patterns)
{
    write_file(input, barcodes, realBarcodes, patterns);

    barcodes.clear();
    realBarcodes.clear();
}

void split_barcodes(input& input, BarcodeMappingVector& barcodes, BarcodeMappingVector& realBarcodes, BarcodePatternVectorPtr barcodePatterns, fastqStats& fastqStatsFinal,
                        const std::vector<std::pair<std::string, char> >& patterns)
{
    //read all fastq lines into str vector
    gzFile fp;
    kseq_t *ks;
    fp = gzopen(input.inFile.c_str(),"r");
    if(NULL == fp){
        fprintf(stderr,"Fail to open file: %s\n", input.inFile.c_str());
    }
    int totalReads = 0;
    unsigned char buffer[1000];
    while(!gzeof(fp))
    {
        gzread(fp, buffer, 999);
        for (const char &c : buffer) 
        {
            if ( c == '\n' )
            {
                ++totalReads;
            }
        }
    }
    totalReads = totalReads/4;
    gzrewind(fp);

    ks = kseq_init(fp);
    std::vector<std::string> fastqLines;
    fastqStats emptyStats = fastqStatsFinal; //an empty statistic with already the right keys for each barcode in the dictionary, used to initialize the stats for each thread

    bool readsLeft = true;
    int totalCurrentReadNum = 0;
    while(readsLeft)
    {
        //push a batch of reads into a temporary vector
        fastqLines.clear();
        int fastqReadThreashold = input.fastqReadBucketSize;
        int numFastqReads = 0;
        while(numFastqReads<fastqReadThreashold)
        {
            if(kseq_read(ks) < 0)
            {
                readsLeft = false;
                numFastqReads = fastqReadThreashold;
                continue;
            }
            fastqLines.push_back(std::string(ks->seq.s));
            ++numFastqReads;
            ++totalCurrentReadNum;

            double perc = totalCurrentReadNum/ (double)totalReads;
            printProgress(perc);
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
        std::vector<fastqStats> statsThreadList(input.threads, emptyStats); // vector of stats to be filled
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

        writeCurrentResults(barcodes, realBarcodes, fastqStatsFinal, input, patterns);
    }
    std::cout << "\nMATCHED: " << fastqStatsFinal.perfectMatches << " | MODERATE MATCH: " << fastqStatsFinal.moderateMatches
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
                if(!(c=='A' || c=='T' || c=='G' || c=='C' ||
                    c=='a' || c=='t' || c=='g' || c=='c' ||
                    c=='N' || c=='X'))
                {
                    std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C,N)\n";
                    if(c==' ' || c=='\t' || c=='\n')
                    {
                        std::cerr << "PARAMETER ERROR: Detected a whitespace in sequence; remove it to continue!\n";
                    }
                    exit(1);
                }
                if(c!='N'){nonConstantSeq=false;}

                //determine pattern type and throw error
                if( (patternType!='c') && (c!='N') && (c!='X') )
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
                if(!(c=='A' || c=='T' || c=='G' || c=='C' ||
                        c=='a' || c=='t' || c=='g' || c=='c'))
                        {
                        std::cerr << "PARAMETER ERROR: a barcode sequence in barcode file is not a base (A,T,G,C)\n";
                        if(c==' ' || c=='\t' || c=='\n')
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
        barcodeFile.close();
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

void initializeOutput(std::string output, const std::vector<std::pair<std::string, char> > patterns)
{
    //remove output
    std::string outputStats;
    std::string outputReal;
    std::string outputMapped;
    std::size_t found = output.find_last_of("/");
    if(found == std::string::npos)
    {
        outputMapped = "BarcodeMapping_" + output;
        outputStats = "StatsBarcodeMappingErrors_" + output;
        outputReal = "RealBarcodeSequence_" + output;
    }
    else
    {
        outputMapped = output.substr(0,found) + "/" + "BarcodeMapping_" + output.substr(found+1);
        outputStats = output.substr(0,found) + "/" + "StatsBarcodeMappingErrors_" + output.substr(found+1);
        outputReal = output.substr(0,found) + "/" + "RealBarcodeSequence_" + output.substr(found+1);
    }
    // remove outputfile if it exists
    std::remove(outputMapped.c_str());
    std::remove(outputStats.c_str());
    std::remove(outputReal.c_str());

    //write header line
    std::ofstream outputFile;   
    outputFile.open (outputMapped, std::ofstream::app);
    for(int i =0; i < patterns.size(); ++i)
    {
        outputFile << patterns.at(i).first;
        if( i!=(patterns.size() - 1) )
        {
           outputFile << "\t";
        }
    }
    outputFile << "\n";
    outputFile.close();
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
        initializeOutput(input.outFile, patterns);
        initializeStats(fastqStats, barcodePatterns);
        split_barcodes(input, mappedBarcodes, realBarcodes, barcodePatterns, fastqStats, patterns);
        writeStats(input.outFile, fastqStats);
    }
 
    return EXIT_SUCCESS;
}