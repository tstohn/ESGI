#include <iostream>
#include <fstream>
#include <string>

//stores all the input parameters
struct input{
    std::string inFile;
    std::string outFile;

    std::string barcodeFile;
    std::string mismatchLine;
    std::string patternLine;

    int threads = 5;
};

struct fastqStats{
    int perfectMatches = 0;
    int noMatches = 0;
    int moderateMatches = 0;
};