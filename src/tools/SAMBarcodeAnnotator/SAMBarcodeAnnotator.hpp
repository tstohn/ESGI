#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

class SAMBarcodeAnnotator {
public:
    // Constructor
    SAMBarcodeAnnotator(const std::string& barcodeFile, const std::string& samFile);

    // Main function to perform annotation
    void annotate();

private:
    std::string barcodeFile;
    std::string samFile;

    // Hash map storing read names -> barcode information (all columns)
    std::unordered_map<std::string, std::vector<std::string>> readToBarcodeMap;

    // Hash map storing read names -> mapped gene
    std::unordered_map<std::string, std::string> readToGeneMap;

    // Methods
    void readBarcodeFile();
    void processSAMFile();
    void writeAnnotatedSAMFile();
    void writeAnnotatedBarcodeFile();
};

