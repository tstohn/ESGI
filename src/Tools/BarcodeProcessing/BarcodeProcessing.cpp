#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <vector>
#include <thread>
#include<pthread.h>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "dataTypes.hpp"
#include "helper.hpp"

using namespace boost::program_options;
pthread_mutex_t lock;

/**
 * HOW TO:
 * - generates two files: 
 *          - firstly collapse all UMI_barcodes and have a file: UMI_idx, Ab_idx, cell_idx
 *          - secondly collapse all UMI_idxs: Ab_idx, Ab_count, cell_idx
 * 
 *  - get all unique UMIs
 *  - cluster UMIs that r similar to the same group of UMIs
 * => make a file where UMI column is the 'expected' UMI, from this file we can calculate the saturation...
 * 
 *  - delete doublicate reads (based on UMI and BC 1-4)
 *  - for the same cell (based on BC 1-4) count the number of UMIs per AB_barcode and write a new column AB_count,
 *    so that we end up with a file Ab_id, sample_id, Ab_count
 * => generate the file format for our ScRNA_seq_Normalization pipeline
 * 
 * */

struct dataLine
{
    char* umi_seq;
    char* ab_seq;
    char* cell_seq;
    
    //dataLine():val(0), i(0), j(0){}

};

class UmiData
{


    public:
        ~UmiData()
        {
            data.clear();
            uniqueChars.~UniqueCharSet();
        }

        void parseFile(const std::string fileName, const int& thread)
        {
            int totalReads = totalNumberOfLines(fileName);
            int currentReads = 0;
            //open gz file
            std::ifstream file(fileName, std::ios_base::in | std::ios_base::binary);
            boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
            inbuf.push(boost::iostreams::gzip_decompressor());
            inbuf.push(file);
            std::istream instream(&inbuf);

            //run multiple threads on bucket of barcode lines
            std::vector<std::thread> workers;
            for (int i = 0; i < thread; ++i) 
            {
                workers.push_back(std::thread(&UmiData::parseBarcodeLines, this, &instream, totalReads, std::ref(currentReads) ));
            }
            for (std::thread &t: workers) 
            {
                if (t.joinable()) {
                    t.join();
                }
            }

            //Cleanup
            file.close();
        }

    private:

        dataLine extraDataLineFromString(const std::string& line)
        {
            dataLine dataTmp;

            //split the line into barcodes
            std::vector<std::string> result;
            std::stringstream ss;
            ss.str(line);
            while( ss.good() )
            {
                std::string substr;
                getline( ss, substr, ',' );
                result.push_back( substr );
            }

            return dataTmp;
        }

        void parseBarcodeLines(std::istream* instream, const int& totalReads, int& currentReads)
        {
            std::string line;
            pthread_mutex_lock(&lock);
            if(!std::getline(*instream, line))
            {
                pthread_mutex_unlock(&lock);
                dataLine datalineTmp = extraDataLineFromString(line);   
                pthread_mutex_lock(&lock);
             
            }

            double perc = currentReads/ (double)totalReads;
            ++currentReads;
            printProgress(perc);
            pthread_mutex_unlock(&lock);



            dataLine dataTmp = extraDataLineFromString(line);
            data.push_back(dataTmp);
        }

        std::vector<dataLine> data;
        UniqueCharSet uniqueChars;

};

class AbData
{
    public:

        AbData(UmiData umiData)
        {
            
        }


        void writeFile(const std::string outFile)
        {

        }
};

bool parse_arguments(char** argv, int argc, std::string& inFile,  std::string& outFile, int& threats)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            ("input,i", value<std::string>(&inFile)->required(), "directory of files or single file in fastq(.gz) format")
            ("output,o", value<std::string>(&outFile)->required(), "output file with all split barcodes")

            ("thread,t", value<int>(&threats)->default_value(5), "number of threads")
            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";
            std::cout << "EXAMPLE CALL:\n ./bin/processing -i <inFile>\n";
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


int main(int argc, char** argv)
{
    std::string inFile;
    std::string outFile;
    int thread;
    parse_arguments(argv, argc, inFile, outFile, thread);
    UmiData data;
    data.parseFile(inFile, thread);

    AbData abData(data);
    abData.writeFile(outFile);

    return(EXIT_SUCCESS);
}