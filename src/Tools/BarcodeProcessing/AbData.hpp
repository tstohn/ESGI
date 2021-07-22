#pragma once

#include "UmiData.hpp"

struct abLine
{
    const char* ab_seq;
    const char* cell_seq;
    int ab_cout;

};
typedef std::shared_ptr<dataLine> dataLinePtr;
 

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
