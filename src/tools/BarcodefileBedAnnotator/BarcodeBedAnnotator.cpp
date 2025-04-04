#include "BarcodeBedAnnotator.hpp"
#include <regex>

// Constructor
BarcodeBedAnnotator::BarcodeBedAnnotator(const std::string &barcodeFile, const std::string &bedFile, const int featureCol)
    : barcodeFile(barcodeFile), bedFile(bedFile), featureCol(featureCol) {}

// Main function to execute annotation
void BarcodeBedAnnotator::annotate() 
{
    //read barcodes into hash map
    readBedFile();
    //process line by line and add a new gene column to <barcodeFile>_annotated.tsv
    processBarcodeFile();
}

std::string BarcodeBedAnnotator::getNthElement(const std::string& line, int n) 
{
    std::stringstream ss(line);
    std::string element;
    int count = 0;

    while (std::getline(ss, element, '\t')) 
    {
        if (count == n) 
        {
            return element;
        }
        count++;
    }

    std::cout << "Could not find " << std::to_string(n) << "th element in line\n";
    exit(EXIT_FAILURE);
}

// Reads the bed file (rna-read annotations) into a hashmap (readName -> feature)
void BarcodeBedAnnotator::readBedFile() 
{
    std::ifstream bedFileStream(bedFile); 
    if (!bedFileStream) 
    {
        std::cerr << "Error: Unable to open BED file " << bedFile << ".\n";
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(bedFileStream, line)) 
    {
        //skip empty and comment lines
        if (line.empty() || line[0] == '#') continue; 

        std::istringstream stream(line);
        std::vector<std::string> columns;
        std::string token;

        // Split line into columns
        while (std::getline(stream, token, '\t')) 
        {
            columns.push_back(token);
        }

        // Ensure at least featureColumn presence
        if (columns.size() < featureCol+1 || columns.size() < 4)
        {
            std::cout << "line in bed file does not have enough columns to extract feature and read name!\n";
            std::cout << "Detecting only " << std::to_string(columns.size()) << " columns. " << std::to_string(featureCol+1) << " expected\n";
            std::cout << line << "\n";
            continue;
        }

        //name is in 4th column in bed file
        std::string readName = columns.at(3);
        // Extract the feature from the target column
        std::string feature = columns.at(featureCol);

        if (!feature.empty()) 
        {
            readnameToFeatureMap.insert(std::make_pair(readName, feature));
        }
    }
}

// Reads the barcode tsv file, annotates it, and extracts gene mappings
void BarcodeBedAnnotator::processBarcodeFile() 
{
    std::ifstream barcodeFileStream(barcodeFile);
    if (!barcodeFileStream) 
    {
        std::cerr << "Error: Unable to open barcode file: " << barcodeFile << std::endl;
        return;
    }
    std::string outputFile;
    size_t pos = barcodeFile.rfind(".tsv");
    if (pos == std::string::npos) 
    {
        std::cerr << "Error: Unable to open barcode file: " << barcodeFile << std::endl;
        exit(EXIT_FAILURE);
    }
    outputFile = barcodeFile.substr(0, pos) + "_annotated" + barcodeFile.substr(pos);
    std::ofstream annotatedFileStream(outputFile);
    if (!annotatedFileStream) 
    {
        std::cerr << "Error: Unable to open output (annotated barcode) file: " << outputFile << std::endl;
        return;
    }

    std::string line;
    std::getline(barcodeFileStream, line);
    annotatedFileStream << line << "\t" << "FEATURE" << "\n";  // Write to output file

    while (std::getline(barcodeFileStream, line)) 
    {
        //read name is in the 1st columns
        std::string readName = getNthElement(line, 0);
        std::string feature = "UNMAPPED";

        auto it = readnameToFeatureMap.find(readName);
        if (it != readnameToFeatureMap.end()) 
        {
            feature = it->second; //get feature name if found
        }

        annotatedFileStream << line << "\t" << feature << "\n";  // Write to output file
    }

    barcodeFileStream.close();
    annotatedFileStream.close();
}
