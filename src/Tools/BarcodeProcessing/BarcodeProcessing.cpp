#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <vector>

#include "dataTypes.hpp"

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
    const char* umi_seq;
    const char* ab_seq;
    const char* cell_seq;
};

class UmiData
{


    public:
        ~UmiData()
        {
            data.clear();
            uniqueChars.~UniqueCharSet();
        }

        void parseFile(const std::string fileName)
        {
            
        }

        void addBarcodeLine(const std::string line)
        {

        }

    private:
        std::vector<dataLine> data;
        UniqueCharSet uniqueChars;

};



int main(int argc, char** argv)
{


    return(EXIT_SUCCESS);
}