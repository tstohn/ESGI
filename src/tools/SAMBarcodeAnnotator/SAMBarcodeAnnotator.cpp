#include "SAMBarcodeAnnotator.hpp"

// Constructor
SAMBarcodeAnnotator::SAMBarcodeAnnotator(const std::string &barcodeFile, const std::string &samFile)
    : barcodeFile(barcodeFile), samFile(samFile) {}

// Main function to execute annotation
void SAMBarcodeAnnotator::annotate() 
{
    //read barcodes into hash map
    readBarcodeFile();
    //process line by line and add barcode information// store genes in hash map
    processSAMFile();

    //write the lines with added barcodes from line above
    writeAnnotatedSAMFile();
    //combine hash map of barcodes with the filled gene hash map
    writeAnnotatedBarcodeFile();
}

// Reads the barcode file into a hashmap (readName -> all barcode columns)
void SAMBarcodeAnnotator::readBarcodeFile() 
{
    std::ifstream barcodeFileStream(barcodeFile);
    if (!barcodeFileStream) {
        std::cerr << "Error: Unable to open barcode file: " << barcodeFile << std::endl;
        return;
    }

    std::string line;
    std::getline(barcodeFileStream, line); // Read header (we will append "Gene" later)

    while (std::getline(barcodeFileStream, line)) 
    {
        std::stringstream ss(line);
        std::string readName;
        std::vector<std::string> barcodeData;
        std::string value;

        ss >> readName; // First column is the read name
        while (ss >> value) 
        {
            barcodeData.push_back(value); // Store all remaining columns
        }

        readToBarcodeMap[readName] = barcodeData;
    }

    barcodeFileStream.close();
}

// Reads the SAM file, annotates it, and extracts gene mappings
void SAMBarcodeAnnotator::processSAMFile() 
{
    std::ifstream samFileStream(samFile);
    if (!samFileStream) {
        std::cerr << "Error: Unable to open SAM file: " << samFile << std::endl;
        return;
    }

    std::vector<std::string> samLines;
    std::string line;

    while (std::getline(samFileStream, line)) 
    {
        // Keep header lines as they are
        if (line[0] == '@') 
        {
            samLines.push_back(line);
            continue;
        }

        std::stringstream ss(line);
        std::string readName, gene;
        ss >> readName; // First field is the read name

        // Extract the gene name (assume it is stored in the third column of SAM)
        std::vector<std::string> fields;
        std::string field;
        while (ss >> field) 
        {
            fields.push_back(field);
        }

        if (fields.size() > 2) 
        { // Ensure there is a gene field (adjust index if necessary)
            gene = fields[2];
        } else 
        {
            gene = "NA"; // Default to NA if no gene found
        }

        // Add barcode annotation to SAM file
        if (readToBarcodeMap.find(readName) != readToBarcodeMap.end()) {
            for (const auto &barcode : readToBarcodeMap[readName]) {
                line += "\t" + barcode;
            }
        } else {
            line += "\tNA"; // No barcode found
        }

        samLines.push_back(line);

        // Store gene mapping for barcode file annotation
        readToGeneMap[readName] = gene;
    }

    samFileStream.close();

    // Write the annotated SAM file
    std::ofstream annotatedSAMFile("annotated_" + samFile);
    if (!annotatedSAMFile) {
        std::cerr << "Error: Unable to open output SAM file for writing" << std::endl;
        return;
    }

    for (const auto &annotatedLine : samLines) {
        annotatedSAMFile << annotatedLine << std::endl;
    }

    annotatedSAMFile.close();
}


void SAMBarcodeAnnotator::writeAnnotatedSAMFile() 
{
}

// Writes the barcode file with added gene information
void SAMBarcodeAnnotator::writeAnnotatedBarcodeFile() {
    std::ifstream barcodeFileStream(barcodeFile);
    std::ofstream annotatedBarcodeFile("annotated_" + barcodeFile);

    if (!barcodeFileStream || !annotatedBarcodeFile) {
        std::cerr << "Error: Unable to open files for annotated barcode output" << std::endl;
        return;
    }

    std::string line;
    std::getline(barcodeFileStream, line);
    annotatedBarcodeFile << line << "\tGene" << std::endl; // Append "Gene" column

    while (std::getline(barcodeFileStream, line)) {
        std::stringstream ss(line);
        std::string readName;
        ss >> readName;

        if (readToGeneMap.find(readName) != readToGeneMap.end()) {
            line += "\t" + readToGeneMap[readName];
        } else {
            line += "\tNA"; // If no gene found, add 'NA'
        }

        annotatedBarcodeFile << line << std::endl;
    }

    barcodeFileStream.close();
    annotatedBarcodeFile.close();
}