#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <zlib.h>
#include <regex>
#include <thread>

#include "Barcode.hpp"
#include "Demultiplexer.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/version.hpp>

/** 
 * @brief A tool to map fastq lines (stitched to one read, use e.g. fastq-join) to a certain barcode pattern,
 * this pattern must be given as an input parameter plus additional files for patterns with multiple
 * possible barcodes
 * 
 * HOW TO ALGORITHM:
    using edlib (bit-parallel algorithm) with maximum MM number and semi-global alignment (allowed deletions in the read-sequence without punishment).
    map pattern-elements one after the other and hand over sub-sequences of length (pattern-element-length + allowed number MM) to
    edlib-align function, then cut the pattern-element at the alignment end before the start of unpunished deletions. The last edit-operations
    are forced to be substitutions in the case where it can be deletions or substitutions. UMIs are just 'cut-out' at the expected length.
    Before aligning we check a hash-map if read-sequence contains a perfect barcode match. Additionally, we end aligning variables barcodes early
    when we fin a match that is below a threshold (e.g., minimum conversion rate/2, bcs. the read-sequence can be in the middle between two
    barcodes, so the conversion rate divided by two is the threshoold where this number of mismatches could define non-unique barcodes).

 * @param <input> input fastq file (gzipped), considers only full length fastq reads, FW/RV must be stitched together in advance
 * @param <output> output extension, that will be added to the output files, see return
 * @param <sequencePattern> a string with all the barcode patterns, each pattern is enclosed by suqare brackets, valid chars are AGTC ofr bases, N for a sequences
 *                          that hold a fixed number of different combinations, must be declared in the barcodeList file, D and X for wildcard sequences, that can
 *                          be anything - where X is used for e.g. UMIs and D can be used for subsequent mapping if 'writeTranscriptomeFastq'-flag (-c true) is defined
 *                          (e.g. for UMIs): e.g.: [XXXXXXXXXX][AGCGTACGCGAGT][NNNNNNNN][AAGCGtAGCTTC][NNNNNNNN] 
 * @param <barcodeList> a file of discrete barcodes, that can be mapped to the [N...] sequences, each row is one pattern, they r in the order as in the sequence pattern and
 *                      barcodes must be comma seperated
 * @param <mismatches> comma seperated list of mismatches, one entry for each sequence pattern declared above: e.g.: 1,2,1,2,1. This also requires a parameter
                       for the 'X'-sequence, although this parameter is at the moment unused, as we have no sequence to alig this one to.
                       Be aware that in paired-end mapping both strands are mapped, and only reads where all overlapping BCs are the same are kept.
                       For the UMI this means the UMI must be 100% preserved, as we do no UMI-mapping step at this point, so the MM for the UMI-position is 
                       completely unused
 *  
 * @return  writes a tab seperated file with barcodes, a file with statistical information, and if flag <-r> set, also a tab seperated file
 *          of corresponding reads with their real sequence, including mismatches.
 *          After running a line of perfect, moderate and mismatches is printed. Each count refers to a whole line (all barcodes matched perfectly,
 *          at least one matched only with mismatches, at least one did not match at all. The number of duplicate barcodes is per barcode, and thus might
 *          be greater than the total number of mismatches.)
 * 
 *          barcode file: Output file starts with 'BarcodeMapping_' and ends with the ending extension declared with the -o parameter
 *          stats file: a list where for each barcode in the fastq file, the number of found mismatches is written, the columns start with zero mismatches
 *          to the maximum number of mismtaches declared in mismatches parameter plus one, this last columns sums up all cases of more than max-num mismatches.
 *          real barcode file: like barcode file with uncorrected sequences, starts with RealBarcodeSequence
 **/

using namespace boost::program_options;

bool parse_arguments(char** argv, int argc, input& input)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input,i", value<std::string>(&(input.inFile))->required(), "single file in fastq(.gz) format or the forward read file, if <-r> is also set for the\
            reverse reads. It is also possible to provide a txt file with fastq-lines only (with no fastq-quality lines. For txt-files only the single-read option\
            with forward read only is supported: -i file.txt)")
            //optional for reverse mapping: no recommended, join reads first
            ("reverse,r", value<std::string>(&(input.reverseFile))->default_value(""), "Use this parameter for paired-end analysis as the reverse read file. <-i> is the forward read in \
            this case.")
            ("independent,d", value<bool>(&(input.independentReverseMapping))->default_value(false),"independent mapping of forward and reverse read. In this case we do not \
            assume the whole pattern is one sequence from 5'->3'. We rather have two seperate reads for FW and RV and we map both reads individually and the reverse\
            read is not a reverse complement of the pattern itself. In this case we must additionally add a read seperator [-] to clarify where FW and RV reads end. \
            Barcodes for the reverse read are then mapped as they are and are not reverse complements of the pattern.")

            ("output,o", value<std::string>(&(input.outPath))->required(), "output directory. All files including failed lines, statistics will be saved here.")
            ("namePrefix,n", value<std::string>(&(input.prefix))->default_value(""), "a prefix for file names. Default uses no prefix.")

            //removed: old parameter for seuqence pattern: its now parsed directly form the demultiplexed output
            /*("sequencePattern,p", value<std::string>(&(input.patternLine))->required(), "pattern for the sequence to match, 
            every substring that should be matched is enclosed with square brackets. N is a barcode match, X is a wild card match 
            and D is a transcriptome read (e.g. cDNA), * is a stop sign (must be enclosed in brackets [*] and then mapping stops at this position 
            on both sides from FW and RV read): [AGCTATCACGTAGC][XXXXXXXXXX][NNNNNN][AGAGCATGCCTTCAG][NNNNNN]")*/

            ("barcodePatternsFile,p", value<std::string>(&(input.barcodePatternsFile))->required(), "patterns for the sequences to match, every substring that should be matched is enclosed with square brackets. \
            Linker sequences of known barcodes can be 'hard-coded', e.g. [ACGTCAG], for variable barcodes one can add a file path, e.g.[data/barcodes.txt] with comma seperated possible barcodes that can be found at this position(the barcodes can be of variable lengths), [10X] is a wild card match with 10 random bases (e.g., for UMIs) and [DNA] is a transcriptome read (e.g. cDNA), [*] is a stop sign/random sequence part that will also not be mapped, and [-] seperates forward and reverse read. ALL signs (also [*] and [-] must be enclosed in brackets). \
            You can supply several rows with various patterns that might all exist in the input fastq.\
            \nSIGN DETAILS: \
            \n[*]: mapping stops at this position on both sides from FW and RV read, completely disregarding any sequence that follows from FW/ RV read. This sign can only be used ONCE in a pattern). E.g.: [AGCTATCACGTAGC][XXXXXXXXXX][BC1.txt][*][AGAGCATGCCTTCAG][BC1.txt]. in the FW read we map only [AGCTATCACGTAGC][XXXXXXXXXX][BC1.txt] and in the reverse read only [AGAGCATGCCTTCAG][BC1.txt].\
            \n[-]: seperates FW and RV reads. Useful if one read contains DNA and (after mapping a barcode in the beginning) the rest of the read should be assign to DNA. E.g.: [BC1.txt][DNA][-][15X][BC1.txt]: extracts BC1 in FW read and assigns the rest of the read to DNA that can be mapped to a reference. \
            \n[DNA]: DNA can only be at the end of a read (we map the FW and RV read from the 5' to the 3' end), this mean that we can have any barcodes before a DNA pattern, but we can not have barcodes after. In other words the DNA pattern must always be at the end of a read and valid patterns must look like this (where ... can be any pattern except [*],[-]): ...[DNA][-][DNA]... \
            \nvalid structures: [BC1.txt][DNA][-][15X][BC1.txt], [BC1.txt][-][DNA][15X][BC1.txt], [DNA][-][15X][BC1.txt], [15X][BC1.txt][-][DNA], [BC1.txt][15X][BC1.txt][DNA] \
            \nNOT valid structures: [BC1.txt][DNA][BC2.txt][-][15X][BC1.txt], [BC1.txt][DNA][15X][BC1.txt]")

            ("mismatchFile,m", value<std::string>(&(input.mismatchFile))->default_value(""), "File with lists of mismatches allowed for each bracket enclosed sequence substring. \
            This should be a comma seperated list of numbers for each substring of the sequence enclosed in squared brackets. E.g.: 2,1,2,1,2. (Also add mismatches for the STOP[*], UMI[X], READSEPERATOR[-] -  \
            this number is not used however, UMIs are aligned in BarcodeProcessing.) We need one line for every line in the barcodePatternsFile.")

            ("threads,t", value<int>(&(input.threads))->default_value(5), "number of threads")
            ("fastqReadBucketSize,s", value<long long int>(&(input.fastqReadBucketSize))->default_value(-1), "number of lines of the fastQ file that should be read into RAM \
            and be processed, before the next fastq read is processed. By default it equal to 100K lines a time.")
            ("writeStats,q", value<bool>(&(input.writeStats))->default_value(false), "writing Statistics about the barcode mapping. This creates three files: \
            ..._Quality_lastPositionMapped.txt stores how often mapping failed at which position for reads that could not be mapped (THIS IS ONLY WRITTEN IF WE HAVE ONLY ONE PATTERN) \
            ..._Quality_typeMM.txt stores for every barcode how often we observed a Subst, Ins, Del \
            ..._Quality_numberMM.txt stores how many mismatches we observed in which barcodes \n")
            ("writeFailedLines,f", value<bool>(&(input.writeFailedLines))->default_value(false), "write failed lines to an extra file.\n")
            ("hamming,H", value<bool>(&(input.hamming))->default_value(false), "Use hamming distance of 1 for all variable barcodes. \
            This way only 1 substitution per barcode is allowed. It is useful for, e.g. 10X data where barcodes can be easily converted and reducing errors \
            to substitutions reduces runtime. \n")

            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";

            std::cout << "###########################################\n";
            std::cout << "EXAMPLE CALL:\n ./bin/demultiplex -i ./src/test/test_data/test_input/testBig.fastq.gz -o ./bin/ -p ./src/test/test_data/test_input/barcodePatternsBig.txt -m ./src/test/test_data/test_input/barcodeMismatchesBig.txt -t 1 -f 1 -q 1 \n";
            std::cout << "For a better understanding of how the barcode and mismatch files should look like, look into ./src/test/test_data/test_input/barcodePatternsBig.txt and ./src/test/test_data/test_input/barcodeMismatchesBig.txt \n";
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

void write_parameter_file(const input& input)
{
    std::string parameterFile = "parameter.ini";
    if(input.prefix !="")
    {
        parameterFile = input.prefix + "_" + parameterFile;
    }

    std::ofstream outFile(input.outPath + "/" + parameterFile);
    if (!outFile) {
        std::cerr << "Could not open parameter file for writing.\n";
        exit(EXIT_FAILURE);
    }

    outFile << "inFile = " << input.inFile << "\n";
    outFile << "outPath = " << input.outPath << "\n";
    outFile << "prefix = " << input.prefix << "\n";
    outFile << "reverseFile = " << input.reverseFile << "\n";
    outFile << "barcodeFile = " << input.barcodeFile << "\n";
    outFile << "patternLine = " << input.patternLine << "\n";
    outFile << "independent reverse read = " << input.independentReverseMapping << "\n";

    outFile << "writeStats = " << (input.writeStats ? "true" : "false") << "\n";
    outFile << "writeFailedLines = " << (input.writeFailedLines ? "true" : "false") << "\n";
    outFile << "writeFilesOnTheFly = " << (input.writeFilesOnTheFly ? "true" : "false") << "\n";

    outFile << "fastqReadBucketSize = " << input.fastqReadBucketSize << "\n";
    outFile << "threads = " << input.threads << "\n";
    
    // Write mismatchFile path and its contents
    outFile << "mismatchFile = " << input.mismatchFile << "\n";
    std::ifstream mismatchIn(input.mismatchFile);
    if (mismatchIn) {
        std::string line;
        while (std::getline(mismatchIn, line)) {
            outFile << "  " << line << "\n";
        }
    } else {
        outFile << "  [Could not read mismatchFile]\n";
    }
    
    // Write patternLine path and its contents
    outFile << "barcodePatternsFile = " << input.barcodePatternsFile << "\n";
    std::ifstream patternIn(input.barcodePatternsFile);
    if (patternIn) {
        std::string line;
        while (std::getline(patternIn, line)) {
            outFile << "  " << line << "\n";
        }
    } else {
        outFile << "  [Could not read barcodePatternsFile file]\n";
    }

    outFile.close();
}

int main(int argc, char** argv)
{

    input input;
    if(parse_arguments(argv, argc, input))
    {
        //check output is a valid directory
        if(! (std::filesystem::exists(input.outPath) && std::filesystem::is_directory(input.outPath)))
        {
            fprintf(stderr,"The output directory (-o) must exist! Please provide a valid directory.\n Fail to find directory: %s\n", input.outPath.c_str());
            exit(EXIT_FAILURE);
        }
        //write parameters to a parameter file
        write_parameter_file(input);

        //set the number of reads in the processing queue by default to 10X number of threads
        if(input.fastqReadBucketSize == -1)
        {
            input.fastqReadBucketSize = input.threads * 100000;
        }

        // run demultiplexing
        if(!input.reverseFile.empty())
        {
            //run in paired-end mode (allowing only fastq(.gz) format)
            if(!(endWith(input.inFile, "fastq") || endWith(input.inFile, "fastq.gz")))
            {
                std::cout << "Wrong file format for forward-read file <-i>!\n";
                exit(EXIT_FAILURE);
            }
            if(!(endWith(input.reverseFile, "fastq") || endWith(input.reverseFile, "fastq.gz")))
            {
                std::cout << "Wrong file format for reverse-read file <-r>!\n";
                exit(EXIT_FAILURE);
            }

            Demultiplexer<MapEachBarcodeSequentiallyPolicyPairwise, ExtractLinesFromFastqFilePolicyPairedEnd> mapping;
            mapping.run(input);
        }
        else if(endWith(input.inFile, "fastq") || endWith(input.inFile, "fastq.gz"))
        {
            Demultiplexer<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy> mapping;
            mapping.run(input);
        }
        else if(endWith(input.inFile, "txt"))
        {
            Demultiplexer<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy> mapping;
            mapping.run(input);
        }
        else
        {
            fprintf(stderr,"Input file must be of format: <.fastq> | <.fastq.gz> | <.txt>!!!\nFail to open file: %s\n", input.inFile.c_str());
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
