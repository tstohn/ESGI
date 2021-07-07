#include <iostream>
#include <string>
#include <zlib.h>
#include <regex>
#include <thread>

#include "seqtk/kseq.h"

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

KSEQ_INIT(gzFile, gzread)

struct barcode{
    std::string antibody; //BARCODE 1
    std::string sample; //BARCODE 2
    std::string umi;
};
typedef std::vector<barcode> barcodeVector;
using namespace boost::program_options;

struct input{
    std::string inFile;
    std::string outFile;

    std::string wellBarcodes;
    std::string abBarcodes;

    std::string anchor;
    int mb = 1;
    int ma = 2;
    int threads = 5;
};

struct fastqStats{
    int perfectMatches = 0;
    int noMatches = 0;
};

bool parse_arguments(char** argv, int argc, input& input)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input,i", value<std::string>(&(input.inFile))->required(), "directory of files or single file in fastq(.gz) format")
            ("output,o", value<std::string>(&(input.outFile))->required(), "output file with all split barcodes")

            ("anchor,a", value<std::string>(&(input.anchor))->required(), "anchor sequence")
            ("well,w", value<std::string>(&(input.wellBarcodes))->required(), "well barcodes")
            ("antibody,ab", value<std::string>(&(input.abBarcodes))->required(), "antibody barcodes")

            ("mb", value<int>(&(input.mb))->default_value(1), "mismatches allowed in barcodes")
            ("ma", value<int>(&(input.ma))->default_value(2), "mismatches allowed in anchor sequence")
            ("t", value<int>(&(input.threads))->default_value(5), "number of threads")

            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";
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

bool split_sequence(const std::string seq, input* input, barcode& barcode, fastqStats& stats)
{
    //match regex of line
    std::smatch sm;
    std::string expressionString = "[ACTGN]+([ACTGN]{15,15})([ACTGN]{10,10})" + input->anchor + "([ACTGN]{10,10})[ACTGN]+";
    std::regex expression(expressionString);

    std::string line = seq;
    std::regex_match(line,sm,expression);

    //check of mismatches
    if(sm.size()<4)
    {
        //minor mismatches that are allowed per mb and ma
        ++stats.noMatches;
        //bigger mismatches
        return false;
    }
    //perfect matches
    else
    {
        assert(sm.size()==4);
        barcode.umi = sm[1];
        barcode.antibody = sm[2];
        barcode.sample = sm[3];
        ++stats.perfectMatches;

        return true;
    }
}

void write_file(std::string output, barcodeVector barcodes)
{

}

void generate_barcodes(std::vector<std::string> fastqLines, input* input, barcodeVector& barcodes, fastqStats& stats)
{
    barcode barcode;
    for(const std::string& line : fastqLines)
    {
        if(split_sequence(line, input, barcode, stats))
        {
            barcodes.push_back(barcode);
        }
    }
}

barcodeVector iterate_over_fastq(input& input)
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
    std::vector<std::vector<std::string> > fastqLinesVector;
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
    std::vector<barcodeVector> barcodesThreadList(input.threads);
    std::vector<fastqStats> statsThreadList(input.threads);
    std::vector<std::thread> workers;
    //final data variables
    barcodeVector barcodeVectorFinal;
    fastqStats fastqStatsFinal;
    for (int i = 0; i < input.threads; ++i) {
        workers.push_back(std::thread(generate_barcodes, fastqLinesVector.at(i), &input, std::ref(barcodesThreadList.at(i)), std::ref(statsThreadList.at(i))));
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
        barcodeVectorFinal.insert(barcodeVectorFinal.end(), barcodesThreadList.at(i).begin(), barcodesThreadList.at(i).end());
        fastqStatsFinal.perfectMatches += statsThreadList.at(i).perfectMatches;
        fastqStatsFinal.noMatches += statsThreadList.at(i).noMatches;
    }
    std::cout << "MATCHED: " << fastqStatsFinal.perfectMatches << " | MISMATCHED: " << fastqStatsFinal.noMatches << "\n";
    kseq_destroy(ks);
    gzclose(fp);
    return barcodeVectorFinal;
}

int main(int argc, char** argv)
{
    input input;
    parse_arguments(argv, argc, input);
    barcodeVector barcodes = iterate_over_fastq(input);
    write_file(input.outFile, barcodes);
 
    return EXIT_SUCCESS;
}