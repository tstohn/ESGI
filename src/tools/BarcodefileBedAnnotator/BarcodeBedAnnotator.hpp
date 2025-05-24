#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

class BarcodeBedAnnotator {
public:
    // Constructor
    BarcodeBedAnnotator(const std::string& barcodeFile, const std::string& bedFile, const int featureCol);

    // Main function to perform annotation
    void annotate();

private:
    std::string barcodeFile;
    std::string bedFile;
    int featureCol;

    std::unordered_map<std::string, std::string> readnameToFeatureMap;

    // Methods
    std::string getNthElement(const std::string& line, int n);
    //read the bed file creating the map of readname to feature
    void readBedFile();

    //processes each demultiplexed barcode file (tsv)
    void processBarcodeFile();
};

