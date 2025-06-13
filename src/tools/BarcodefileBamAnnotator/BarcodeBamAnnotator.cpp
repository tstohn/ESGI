#include "BarcodeBamAnnotator.hpp"

#include <htslib/sam.h>
#include <regex>

// Constructor
BarcodeBamAnnotator::BarcodeBamAnnotator(const std::string &barcodeFile, const char* bamFile, const char* featureTag)
    : barcodeFile(barcodeFile), bamFile(bamFile), featureTag(featureTag) {}

// Main function to execute annotation
void BarcodeBamAnnotator::annotate() 
{
    //read barcodes into hash map
    readBamFile();
    //process line by line and add a new gene column to <barcodeFile>_annotated.tsv
    processBarcodeFile();
}

std::string BarcodeBamAnnotator::getNthElement(const std::string& line, int n) 
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
void BarcodeBamAnnotator::readBamFile() 
{
    // Open BAM file
    samFile* samFilePtr = sam_open(bamFile, "r");
    if(!samFilePtr) 
    {
        std::cerr << "Failed to open BAM file\n";
        exit(EXIT_FAILURE);
    }

    // Read header
    bam_hdr_t* header = sam_hdr_read(samFilePtr);
    if (!header) 
    {
        std::cerr << "Failed to read BAM header\n";
        sam_close(samFilePtr);
        exit(EXIT_FAILURE);
    }

    // Allocate alignment
    bam1_t* aln = bam_init1();

    // Read alignments
    while(sam_read1(samFilePtr, header, aln) >= 0) 
    {
        // Read name
        std::string readName = bam_get_qname(aln);

        // Get feature-tag, per default the GX tag
        uint8_t* featureTagEntry = bam_aux_get(aln, featureTag);
        std::string feature = "";
        if(featureTagEntry && *featureTagEntry == 'Z') 
        {
            feature = bam_aux2Z(featureTagEntry);
        }

        // Output
        readnameToFeatureMap.insert(std::make_pair(readName, feature));
    }

    // Cleanup
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(samFilePtr);
}

// Reads the barcode tsv file, annotates it, and extracts gene mappings
void BarcodeBamAnnotator::processBarcodeFile() 
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
        std::string feature;

        auto it = readnameToFeatureMap.find(readName);
        if ( (it != readnameToFeatureMap.end()) && (it->second)!= "-" && (it->second)!= "")
        {
            annotatedFileStream << line << "\t" << it->second << "\n";  // Write to output file
        }
    }

    barcodeFileStream.close();
    annotatedFileStream.close();
}
