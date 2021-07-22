#pragma once

#include "UmiData.hpp"

struct abLine
{
    const char* ab_seq;
    const char* cell_seq;
    int ab_cout = 0;

}; 

class AbData
{
    public:

        AbData()
        {
            uniqueChars = std::make_shared<UniqueCharSet>();
        }

        inline void setUniqueCharSet(std::shared_ptr<UniqueCharSet> newUniqueChars)
        {
            uniqueChars = newUniqueChars;
        }

        inline void addLine(const abLine& newAbLine)
        {
            abLineVector.push_back(newAbLine);
        }

        inline std::vector<abLine> getData() const
        {
            return abLineVector;
        }

    private:
        std::shared_ptr<UniqueCharSet> uniqueChars;
        std::vector<abLine> abLineVector;
};
