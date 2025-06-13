#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

class BarcodeBamAnnotator {
public:
    // Constructor
    BarcodeBamAnnotator(const std::string& barcodeFile, const char* bamFile, const char* featureTag);

    // Main function to perform annotation
    void annotate();

private:
    std::string barcodeFile;
    const char* bamFile;
    const char* featureTag;

    std::unordered_map<std::string, std::string> readnameToFeatureMap;

    // Methods
    std::string getNthElement(const std::string& line, int n);
    //read the bam file creating the map of readname to feature
    void readBamFile();

    //processes each demultiplexed barcode file (tsv)
    void processBarcodeFile();
};

