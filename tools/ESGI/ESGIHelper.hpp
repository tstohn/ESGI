// THIS FILE CONTAINS ALL THE TOOL-CALL NEEDED FOR ESGI:
// like demultiplexing, counting, etc
#pragma once

#include <cstdlib>
#include <string>
#include <iostream>
#include <filesystem>

#include "Barcode.hpp"
#include "Demultiplexer.hpp"
#include "BarcodeProcessingHandler.hpp"
#include <BarcodeMapping.hpp>

//check if we have patterns like -,*. These are not included in the demultiplexed output and therefore
//indices for counting must be adjusted. Additioanly, these two patterns can ONLY BE at the end of a sequence, 
//e.g., we can t have an arbitrary pattern * at ther beginning of a read
inline int get_special_pattern_pos(const std::vector<std::string>& patterns)
{
    int posDash = -1;
    int posStar = -1;

    // find positions, ensure the patterns are only used once
    for (size_t i = 0; i < patterns.size(); ++i) 
    {
        if (patterns[i] == "-" && posDash==-1) {posDash = static_cast<int>(i);}
        else if (patterns[i] == "[-]" && posDash!=-1){std::cerr << "The pattern [-] can only be used ONCE in a pattern! (It denotes the end of fw/rv reads, e.g., for cases where there is no overlap between reads\n";}
        if (patterns[i] == "*" && posStar==-1) {posStar = static_cast<int>(i);}
        else if (patterns[i] == "[*]" && posDash!=-1){std::cerr << "The pattern [*] can only be used ONCE in a pattern! Everything 'behind' this pattern in reading direction is not considered anymore for demultiplexing\n";}
    }

    // Case 1: neither found
    if (posDash == -1 && posStar == -1) return-1;

    // Case 2: only one found
    if (posDash != -1) return posDash;
    if (posStar != -1) return posStar;

    // Case 3: both found - throw error
    if (posStar != -1 && posDash != -1)
    {
        std::cerr << "Error: [-] and [*] can not be used together. You can only use one of them! [*] already implicitely denotes the end of the strand, as nothing is demultiplexed beyond this sign\
        making an additional [-] obsolete\n";
    }

    return -1;
}
//get the position of the DNA pattern
inline int get_DNA_pattern_pos(const std::vector<std::string>& patterns)
{
    auto it = std::find(patterns.begin(), patterns.end(), "DNA");
    return( (it == patterns.end() ? -1 : static_cast<int>(it - patterns.begin())) );
}

//obsolete function, still present for potential future usage, but ESGI DOES NOT call demultiplex but 
//has it integrated and now does run_demultiplex
inline bool call_demultiplex(const ESGIConfig& config)
{
    std::filesystem::path toolPath = std::filesystem::path("bin") / "demultiplex";

    // Build command string with parameters
    std::string cmd = toolPath.string() +
        " --input " + config.forward +
        " --output " + config.output +
        " --threads " + std::to_string(config.threads) +
        " --barcodePatternsFile " + config.patternFile +
        " --mismatchFile " + config.mismatchesFile +
        " --independent " + std::to_string(config.independent) +
        " --hamming " + std::to_string(config.hamming) + 
        " --writeFailedLines " + std::to_string(config.writeFailedLines) + 
        " --writeStats " + std::to_string(config.writeStats);

    // optional arguments
    if(config.reverse.has_value()){cmd += " --reverse " + config.reverse.value();}
    if(config.prefix.has_value()){cmd += " --namePrefix " + config.prefix.value();}
    if(config.fastqReadBucketSize.has_value()){cmd += " --fastqReadBucketSize " + std::to_string(config.fastqReadBucketSize.value());}

    std::cout << "Running demultiplexing as: " << cmd.c_str() << "\n";
    int ret = std::system(cmd.c_str());
    if (ret != 0) 
    {
        std::cout << "DEMULTIPLEX failed with code " << ret << "\n";
        return false;
    }

    return true;
}

inline bool run_demultiplex(const ESGIConfig& config)
{
    //set the arguments of demultiplex from ESGIConfig - what otherwise the parse_arguments function did with the input parameters
    input input;

    input.inFile = config.forward;
    input.independentReverseMapping = config.independent;
    input.outPath = config.output;
    input.barcodePatternsFile =  config.patternFile;
    input.mismatchFile = config.mismatchesFile;
    input.threads = config.threads;
    input.hamming = config.hamming;
    input.writeFailedLines = config.writeFailedLines;
    input.writeStats = config.writeStats;

    // optional arguments
    if(config.reverse.has_value()){input.reverseFile = config.reverse.value();}
    if(config.prefix.has_value()){input.prefix = config.prefix.value();}
    if(config.fastqReadBucketSize.has_value()){input.fastqReadBucketSize = config.fastqReadBucketSize.value();}

    //BELOW CODE IS ESSENTIALLY THE DEMULTIPLEXING MAIN
    try
    {
        //check output is a valid directory
        if(! (std::filesystem::exists(input.outPath) && std::filesystem::is_directory(input.outPath)))
        {
            fprintf(stderr,"The output directory (-o) must exist! Please provide a valid directory.\n Fail to find directory: %s\n", input.outPath.c_str());
            exit(EXIT_FAILURE);
        }

        //set the number of reads in the processing queue by default to 10X number of threads
        if(input.fastqReadBucketSize == -1)
        {
            input.fastqReadBucketSize = input.threads * 100000;
        }

        // run demultiplexing
        if(!input.reverseFile.empty())
        {
            //run in paired-end mode (allowing only fastq(.gz) format)
            if(!(endWith(input.inFile, "fastq") || endWith(input.inFile, "fastq.gz")))
            {
                std::cout << "Wrong file format for forward-read file <-i>!\n";
                exit(EXIT_FAILURE);
            }
            if(!(endWith(input.reverseFile, "fastq") || endWith(input.reverseFile, "fastq.gz")))
            {
                std::cout << "Wrong file format for reverse-read file <-r>!\n";
                exit(EXIT_FAILURE);
            }

            Demultiplexer<MapEachBarcodeSequentiallyPolicyPairwise, ExtractLinesFromFastqFilePolicyPairedEnd> mapping;
            mapping.run(input);
        }
        else if(endWith(input.inFile, "fastq") || endWith(input.inFile, "fastq.gz"))
        {
            Demultiplexer<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromFastqFilePolicy> mapping;
            mapping.run(input);
        }
        else if(endWith(input.inFile, "txt"))
        {
            Demultiplexer<MapEachBarcodeSequentiallyPolicy, ExtractLinesFromTxtFilesPolicy> mapping;
            mapping.run(input);
        }
        else
        {
            fprintf(stderr,"Input file must be of format: <.fastq> | <.fastq.gz> | <.txt>!!!\nFail to open file: %s\n", input.inFile.c_str());
            exit(EXIT_FAILURE);
        }
    }
    catch (const std::exception& e) 
    {
        std::cerr << "Demultiplexing failed with error: " << e.what() << "\n";
        return false;
    }

    return true;
}

// --- capture helper (kept minimal) ---
static std::string run_system_capture(const std::string& cmd) {
    std::string tmp = "star_ver_" + std::to_string(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()
    ) + ".txt";
    std::string full = cmd + " > \"" + tmp + "\" 2>&1";
    if (std::system(full.c_str()) == -1) return std::string();
    std::ifstream f(tmp, std::ios::in | std::ios::binary);
    std::string s((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    f.close();
    std::remove(tmp.c_str());
    return s;
}

// parse "A.B.C[letter]" from an arbitrary string (finds first number run)
static bool parse_version_triplet(const std::string& src,
                                  int& A, int& B, int& C, int& S /*-1 no suffix*/) {
    A = B = C = 0; S = -1;
    size_t i = src.find_first_of("0123456789");
    if (i == std::string::npos) return false;

    auto read_int = [&](size_t& p, int& out)->bool {
        if (p >= src.size() || !std::isdigit((unsigned char)src[p])) return false;
        long v = 0;
        while (p < src.size() && std::isdigit((unsigned char)src[p])) {
            v = v*10 + (src[p]-'0'); ++p;
        }
        out = (int)v; return true;
    };

    if (!read_int(i, A)) return false;
    if (i >= src.size() || src[i++] != '.') return false;
    if (!read_int(i, B)) return false;
    if (i >= src.size() || src[i++] != '.') return false;
    if (!read_int(i, C)) return false;

    if (i < src.size() && std::isalpha((unsigned char)src[i])) {
        S = std::tolower((unsigned char)src[i]) - 'a'; // a=0, b=1, ...
    } else {
        S = -1; // no suffix = newest among same A.B.C
    }
    return true;
}

// return positive if (A1,B1,C1,S1) > (A2,B2,C2,S2), 0 if equal, negative if less
static int compare_versions(int A1,int B1,int C1,int S1, int A2,int B2,int C2,int S2) {
    if (A1 != A2) return (A1 > A2) ? 1 : -1;
    if (B1 != B2) return (B1 > B2) ? 1 : -1;
    if (C1 != C2) return (C1 > C2) ? 1 : -1;
    if (S1 == S2) return 0;
    if (S1 == -1) return 1;        // no suffix newer than any suffix
    if (S2 == -1) return -1;
    return (S1 > S2) ? 1 : -1;     // a(0) < b(1) < c(2) ...
}

// returns 0 if STAR exists and version >= minimumVersion (e.g. "2.7.10b"), else 1
static int check_star_version(const std::string& minimumVersion, std::string& currentVersion) 
{
    std::string out = run_system_capture("STAR --version");
    currentVersion = out;
    if (out.empty()) { std::cerr << "STAR not found or no output\n"; return 1; }

    int a,b,c,s, ra,rb,rc, rs;
    if (!parse_version_triplet(out, a,b,c,s)) {
        std::cerr << "Could not parse STAR version from output: " << out << "\n";
        return 1;
    }
    if (!parse_version_triplet(minimumVersion, ra,rb,rc, rs)) {
        std::cerr << "Invalid required version string: " << minimumVersion << "\n";
        return 1;
    }

    int cmp = compare_versions(a,b,c,s, ra,rb,rc,rs);
    return (cmp >= 0) ? 1 : 0;
}

inline bool call_star(const ESGIConfig& config, ESGIIntermediateFiles& intermediateFiles)
{
    std::string starExecutable = "STAR";
    if(config.STAR.has_value()){starExecutable = config.STAR.value();}

    std::filesystem::path outfilePrefixPath = std::filesystem::path(config.output) / "STAR_";
    std::string outfilePrefix = outfilePrefixPath.string();

    // we need a STAR version > 2.7.10b to also get GX/GN tags
    std::string miniumStar = "2.7.10b";
    std::string currentVersion;
    if(!check_star_version(miniumStar, currentVersion))
    {
        std::cerr << "STAR VERSION is not compatible with ESGI, since we require GX/GN annotations from STAR.\n"
         << " We need a minimum version of STAR " << miniumStar << "\n"
         << " Current version is: " << currentVersion << "\n";
        std::cerr << "Please upgrade STAR to a newer version, or set the path for STAR in the input.ini file with <STAR=...>\n";
        exit(EXIT_FAILURE);
    }

    if(!config.genomeDir.has_value())
    {
        std::cerr << "No parameter <genomeDir> provided, please provide this parameter in the.ini file to enable STAR\n";
    }
    const std::string starCmd = starExecutable +
        " --runThreadN " + std::to_string(config.threads) +
        " --genomeDir " + config.genomeDir.value() +
        " --readFilesIn " + intermediateFiles.demultiplexingOutput + ".fastq" +
        " --outFileNamePrefix " + outfilePrefix +
        " --outSAMtype BAM Unsorted " +
        " --outSAMattributes NH HI AS nM GX GN " +
        " --quantMode TranscriptomeSAM " +
        " --outFilterMultimapNmax 50 " +
        " --outSAMmultNmax 1 --outSAMunmapped Within " +
        " --limitOutSJcollapsed 2000000 " +
        " --twopassMode Basic";

    int rc = std::system(starCmd.c_str());
    if (rc != 0) 
    {
        std::cerr << "STAR failed with code " << rc << "\n";
    }

    return true;
}

inline bool call_annotate(const ESGIConfig& config, ESGIIntermediateFiles& intermediateFiles)
{
    std::string demultiplexingPatternOutput = intermediateFiles.demultiplexingOutput + ".tsv";
    // we hard-coded and gave the STAR-output the prefix <STAR>
    std::string starAlignmentOutput = config.output + "/" + "STAR_Aligned.out.bam";

    const std::string annotateCmd = std::string("./bin/annotate ") +
        " -i " + demultiplexingPatternOutput +
        " -b " + starAlignmentOutput +
        " -f " + config.starFeature;

    int rc = std::system(annotateCmd.c_str());
    if (rc != 0) 
    {
        std::cerr << "annotate failed with code " << rc << "\n";
    }

    return true;
}

inline bool run_rna_mapping(const ESGIConfig& config, ESGIIntermediateFiles& intermediateFiles)
{
    //call STAR
    if(!call_star(config, intermediateFiles))
    {
        std::cerr << "Error: STAR failed\n";
    }

    //call annotate
    if(!call_annotate(config, intermediateFiles))
    {
        std::cerr << "Error: Annotation of demultiplexed reads with STAR-mapping failed\n";
    }

    return true;
}

std::string adjust_position_due_to_special_patterns(const std::string& indexList, const ESGIIntermediateFiles& intermediateFiles)
{
    std::stringstream ss(indexList);
    std::string index;
    std::vector<int> values;

    // Split by commas
    while (std::getline(ss, index, ',')) 
    {
        values.push_back(std::stoi(index));
    }

    // Modify values
    for (int &v : values) 
    {
        //check if we had a DNA pattern that got 'cut out' and therefore need to adjust pattern positions
        int dnaShift = 0;
        if(intermediateFiles.dnaPatternPos >= 0 && v > intermediateFiles.dnaPatternPos){dnaShift = -1;}
        //if the index is behind a -,* we need to substract one bcs this poattern is not in the demultiplex output
        if (v > intermediateFiles.specialPatternPos) v -= 1;
        //also we need to add 1 to every column, bcs we will ahve the additional READNAME column
        v += 1;
        //and finally add the dna shift
        v += dnaShift;
    }

    // Rebuild the string
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) 
    {
        if(i) out << ",";
        out << values[i];
    }

    return(out.str());
}

//old function to call count tool, now integrated into ESGI calling the actual code of count main and performing run_count
//this code is only preserved for potential future usage in other pipelines
inline bool call_count(ESGIConfig& config, const ESGIIntermediateFiles& intermediateFiles, const bool dnaPatternPresent)
{
    std::filesystem::path toolPath = std::filesystem::path("bin") / "count";
    std::filesystem::path patternPath(config.patternFile);
    std::filesystem::path dir = patternPath.parent_path();

    //index re-assignment due to -,*
    if(intermediateFiles.specialPatternPos!= -1 || dnaPatternPresent)
    {
        config.SC_ID = adjust_position_due_to_special_patterns(config.SC_ID, intermediateFiles);
        if(intermediateFiles.dnaPatternPos >=0)
        {
            // if we have a DNA pattern, we firstly set the FEATURE column index to patternlength
            // as the DNA is added as a last column, o-indexed it is at position patternlength
            config.FEATURE_ID = std::to_string(intermediateFiles.patternLength);            
            //however, if there were special patterns (like -,*) that are not in the demultiplexed output we need to adust the FEATURE column index
            //basically: -1 for removing the DNA pattern, -1 if there is a -,* and +1 for READNAME col
            config.FEATURE_ID = adjust_position_due_to_special_patterns(config.FEATURE_ID, intermediateFiles);
        }
        else
        {
            config.FEATURE_ID = adjust_position_due_to_special_patterns(config.FEATURE_ID, intermediateFiles);
        }
        if(config.UMI_ID.has_value()){config.UMI_ID = adjust_position_due_to_special_patterns(config.UMI_ID.value(), intermediateFiles);}
        if(config.ANNOTATION_IDs.has_value()){config.ANNOTATION_IDs = adjust_position_due_to_special_patterns(config.ANNOTATION_IDs.value(), intermediateFiles);}
    }

    // Build command string with parameters
    std::string cmd = toolPath.string() +
        " --input " + intermediateFiles.countingInput +
        " --output " + intermediateFiles.countingOutput +
        " --threads " + std::to_string(config.threads) +
        " --barcodeDir " + dir.string() +
        " --featureIndex " + config.FEATURE_ID +
        " --singleCellIndices " + config.SC_ID + 
        " --umiThreshold " +  std::to_string(config.umiThreshold) +
        " --umiAbundanceThreshold " +  std::to_string(config.umiAbundance) +
        " --umiRemoval " +  std::to_string(config.umiCollapsing) +
        " --hamming " + std::to_string(config.hamming);

    // optional arguments
    if(config.FEATURE_NAMES.has_value()){cmd += " --featureNames " + config.FEATURE_NAMES.value();}
    if(config.ANNOTATION_IDs.has_value()){cmd += " --annotationIdxs " + config.ANNOTATION_IDs.value();}
    if(config.ANNOTATION_NAMES.has_value()){cmd += " --annotationFiles " + config.ANNOTATION_NAMES.value();}
    if(config.UMI_ID.has_value()){cmd += " --umiIndex " + config.UMI_ID.value();}
    if(config.barcodeSharing.has_value()){cmd += " --shareBarcodes " + config.barcodeSharing.value();}
    //if we did overwrite the umiMismatches for intermediate data
    if(intermediateFiles.umiMismatches != -1){cmd += " --mismatches " + std::to_string(intermediateFiles.umiMismatches);}

    std::cout << "Running count as: " << cmd.c_str() << "\n";
    int ret = std::system(cmd.c_str());
    if (ret != 0) 
    {
        std::cout << "COUNT failed with code " << ret << "\n";
        return false;
    }

    return true;
}

std::vector<std::string> split_file_lists(const std::string& s)
{
    std::vector<std::string> result;
    std::regex re("[,\\s]+");  // one or more commas or whitespace

    std::sregex_token_iterator it(s.begin(), s.end(), re, -1);
    std::sregex_token_iterator end;

    for (; it != end; ++it) 
    {
        if (!it->str().empty())
            result.push_back(it->str());
    }
    return result;
}


// generate a dictionary to map sequences to AB(proteins)
std::unordered_map<std::string, std::string > generateProteinDict(std::string abFile, 
                                                                  const std::vector<std::string>& abBarcodes)
{
    std::unordered_map<std::string, std::string > map;
    std::vector<std::string> proteinNames;

    std::ifstream abFileStream(abFile);
    if (!abFileStream.is_open()) 
    {
        std::cerr << "Error: Failed to open antibody file" << abFile << ". Please double check if the file exists.\n";
        exit(EXIT_FAILURE);
    }

    for(std::string line; std::getline(abFileStream, line);)
    {
        std::string delimiter = ",";
        std::string seq;
        size_t pos = 0;
        std::vector<std::string> seqVector;
        while ((pos = line.find(delimiter)) != std::string::npos) 
        {
            seq = line.substr(0, pos);
            line.erase(0, pos + 1);
            proteinNames.push_back(seq);
        }
        seq = line;
        proteinNames.push_back(seq);
    }

    assert(abBarcodes.size() == proteinNames.size());
    for(size_t i = 0; i < abBarcodes.size(); ++i)
    {
        map.insert(std::make_pair(abBarcodes.at(i), proteinNames.at(i)));
    }
    abFileStream.close();

    return map;
}

// generate a dictionary to map sequences to treatments
std::unordered_map<int, std::unordered_map<std::string, std::string >> generateAnnotationDict(const std::vector<std::string>& annotationFiles,
                                                                                              const std::vector<int>& annotationIdxs,
                                                                                              const std::vector< std::vector<std::string>>& annotationBarcodesVector)
{
    std::unordered_map<int, std::unordered_map<std::string, std::string >> allAnnotationsMap;

    if(annotationFiles.size() != annotationIdxs.size())
    {
        std::cerr << "Not every annotation-file has a annotation column or vice versa. There are " << std::to_string(annotationFiles.size()) \
        << " annotation-files with " << std::to_string(annotationIdxs.size()) << " annotation-column. Please provide for every file with annotation information" \
        << " one column to match barcodes.";
        exit(EXIT_FAILURE);
    }

    for(size_t aIdx = 0; aIdx < annotationFiles.size(); ++aIdx)
    {
        std::string aFile = annotationFiles.at(aIdx);
        std::unordered_map<std::string, std::string> barcodeToAnnotationMap;

        std::vector<std::string> annotationNames;
        std::ifstream aFileStream(aFile);
        if (!aFileStream.is_open()) 
        {
            std::cerr << "Error: Failed to open cell-grouping file" << aFile << ". Please double check if the file exists.\n";
            exit(EXIT_FAILURE);
        }

        for(std::string line; std::getline(aFileStream, line);)
        {
            std::string delimiter = ",";
            std::string seq;
            size_t pos = 0;
            std::vector<std::string> seqVector;
            while ((pos = line.find(delimiter)) != std::string::npos) 
            {
                seq = line.substr(0, pos);
                line.erase(0, pos + 1);
                annotationNames.push_back(seq);

            }
            seq = line;
            annotationNames.push_back(seq);
        }
        std::vector<std::string> annotationBarcodes = annotationBarcodesVector.at(aIdx);
        if(annotationNames.size() != annotationBarcodes.size())
        {
            std::cout << "The list of condition-barcodes and condition-names has not the same length!\n";
            std::cout << "Please make sure both lists are of the same size and every condition-barcode is assigned a condition-name\n";
            std::cout << "For annotation file <" << aFile << "> we have " << annotationNames.size() << " annotations and " << annotationBarcodes.size() << " possible barcodes\n";
            exit(EXIT_FAILURE);
        }

        //create mapping of barcode to annotation
        for(size_t i = 0; i < annotationBarcodes.size(); ++i)
        {
            barcodeToAnnotationMap.insert(std::make_pair(annotationBarcodes.at(i), annotationNames.at(i)));
        }
        aFileStream.close();
        //make final col-index to annotation map
        allAnnotationsMap.insert(std::make_pair(annotationIdxs.at(aIdx), barcodeToAnnotationMap) );
    }


    return allAnnotationsMap;
}


inline bool run_count(ESGIConfig& config, const ESGIIntermediateFiles& intermediateFiles, const bool dnaPatternPresent)
{
    std::filesystem::path toolPath = std::filesystem::path("bin") / "count";
    std::filesystem::path patternPath(config.patternFile);
    std::filesystem::path dir = patternPath.parent_path();

    //index re-assignment due to -,*
    if(intermediateFiles.specialPatternPos!= -1 || dnaPatternPresent)
    {
        config.SC_ID = adjust_position_due_to_special_patterns(config.SC_ID, intermediateFiles);
        if(intermediateFiles.dnaPatternPos >=0)
        {
            // if we have a DNA pattern, we firstly set the FEATURE column index to patternlength
            // as the DNA is added as a last column, o-indexed it is at position patternlength
            config.FEATURE_ID = std::to_string(intermediateFiles.patternLength);            
            //however, if there were special patterns (like -,*) that are not in the demultiplexed output we need to adust the FEATURE column index
            //basically: -1 for removing the DNA pattern, -1 if there is a -,* and +1 for READNAME col
            config.FEATURE_ID = adjust_position_due_to_special_patterns(config.FEATURE_ID, intermediateFiles);
        }
        else
        {
            config.FEATURE_ID = adjust_position_due_to_special_patterns(config.FEATURE_ID, intermediateFiles);
        }
        if(config.UMI_ID.has_value()){config.UMI_ID = adjust_position_due_to_special_patterns(config.UMI_ID.value(), intermediateFiles);}
        if(config.ANNOTATION_IDs.has_value()){config.ANNOTATION_IDs = adjust_position_due_to_special_patterns(config.ANNOTATION_IDs.value(), intermediateFiles);}
    }

    //below running the code of count main
    try
    {
        std::string inFile;
        std::string outFile;
        std::string barcodeDir;
        std::string barcodeIndices;
        int thread;
        int umiMismatches = 1;
        double umiThreshold = -1;
        bool umiRemoval = true;
        bool scIdAsString = false;
        std::string fuseBarcodesFile = "";

        //data for protein(ab) and treatment information
        std::string abFile; 
        int featureIdx;
        
        std::vector<std::string> annotationFiles;
        std::vector<int> annotationIdxs;
        std::vector<std::string> abBarcodes;
        //vector of barcode vectors, equivalent to the content in annotationFiles, but that one did not parse the annotations yet
        std::vector<std::vector<std::string>> annotationBarcodesVector;
        std::string umiIdx = "";
        bool hamming = false;
        double umiAbundanceThreshold;

        //set the arguments as in parse_arguments on count main here:
        inFile = intermediateFiles.countingInput;
        outFile =  intermediateFiles.countingOutput;
        thread = config.threads;
        barcodeDir = dir.string();
        barcodeIndices = config.SC_ID;
        try 
        {
            featureIdx = std::stoi(config.FEATURE_ID);
        }
        catch (const std::exception& e) 
        {
            std::cerr << "Invalid Feature ID: " << config.FEATURE_ID << "\n";
        }
        //default value in config is 0.0 (umiThreshold is initially set to -1 - also in count - to distinguish between umiThreshold beeing set in parse_argument or not)
        umiThreshold = config.umiThreshold;
        umiRemoval = config.umiCollapsing;
        scIdAsString = config.SC_ID_string;
        hamming = config.hamming;
        umiAbundanceThreshold = config.umiAbundance;

        if(config.UMI_ID.has_value()){umiIdx = config.UMI_ID.value();}
        //umiMismatches is by default -1 in ESGI, however, in count default is 1 for 1MM, therefore we set it to 1 above and re-set here
        //if umiMismatches != -1
        if(intermediateFiles.umiMismatches != -1){umiMismatches = intermediateFiles.umiMismatches;}
        if(config.FEATURE_NAMES.has_value()){abFile = config.FEATURE_NAMES.value();}
        if(config.ANNOTATION_NAMES.has_value()){annotationFiles = split_file_lists(config.ANNOTATION_NAMES.value());}
        if(config.ANNOTATION_IDs.has_value())
        {
            std::vector<std::string> tmpStringAnnotationIdxs = split_file_lists(config.ANNOTATION_IDs.value());
            std::vector<int> intVec;
            for (const auto& s : tmpStringAnnotationIdxs) 
            {
                try 
                {
                    intVec.push_back(std::stoi(s));
                }
                catch (...) 
                {
                    std::cerr << "Invalid integer: " << s << "\n";
                }
            }
            annotationIdxs = intVec;
        }
        if(config.barcodeSharing.has_value()){fuseBarcodesFile = config.barcodeSharing.value();}
        
        //generate the dictionary of barcode alternatives to idx
        BarcodeInformation barcodeIdData;
        barcodeIdData.hamming = hamming;
        barcodeIdData.umiAbundanceThreshold = umiAbundanceThreshold;

        //get the first line of headers from input file
        std::ifstream file;
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        bool gz = isGzipped(inFile);
        std::istream* instream = openFile(inFile, file, inbuf, gz);
        if (!instream)
        {
            std::cerr << "Error reading input file or file is empty! Please double check if the file exists:" << inFile << std::endl;
            exit(EXIT_FAILURE);
        };

        std::string firstLine;
        if (std::getline(*instream, firstLine)) 
        {  
            bool parseAbBarcodes = true;
            if(abFile.empty()){parseAbBarcodes = false;}
            generateBarcodeDicts(firstLine, barcodeDir, barcodeIndices, barcodeIdData, abBarcodes, parseAbBarcodes, featureIdx, umiRemoval, &annotationBarcodesVector, &annotationIdxs, umiIdx, umiMismatches);
        } 
        else
        {
            std::cerr << "Error reading header line of input file:" << inFile << std::endl;
            exit(EXIT_FAILURE);
        }
        //clean data if necessary
        if (instream != &file) delete instream;
        instream = nullptr;

        BarcodeProcessingHandler dataParser(barcodeIdData);
        //if we have barcodes that must be fused (assign certain barcodes to others, bcs they come, e.g., from the same cell)
        if(fuseBarcodesFile != "")
        {dataParser.parse_barcode_sharing_file(fuseBarcodesFile);}
        //default value in barcodeProcessingHandler is 0.0, no need to set if umiThreshold is -1 here
        if(umiThreshold != -1){dataParser.setUmiFilterThreshold(umiThreshold);}
        dataParser.setumiRemoval(umiRemoval);
        dataParser.setSingleCellIdStyle(scIdAsString);

        //generate dictionaries to map sequences to the real names of Protein/ treatment/ etc...
        std::unordered_map<std::string, std::string > featureMap;
        if(!abFile.empty())
        {
            featureMap = generateProteinDict(abFile, abBarcodes);
        }
        //featureMap is empty if the feature name should stay as they are (no mapping of e.g. barcodes to proteins)
        dataParser.addProteinData(featureMap);
        if(!annotationFiles.empty() && !annotationIdxs.empty())
        {
            //MAPPING ANNOTATION_IDX => (ANNOTATION_BARCODE => ANNOTATION_STRING)
            std::unordered_map< int, std::unordered_map<std::string, std::string >> annotationMap = generateAnnotationDict(annotationFiles, annotationIdxs, annotationBarcodesVector);
            dataParser.addAnnotationData(annotationMap);
        }

        //parse the file of demultiplexed barcodes and
        //add all the data to the Unprocessed Demultiplexed Data (stored in rawData)
        // (AB, annotations already are mapped to their real names, scID is a concatenation of numbers for each barcode in
        //each abrcoding round, seperated by a dot)
        dataParser.parse_barcode_file(inFile);

        //further process the data (correct UMIs, collapse same UMIs, etc.)
        dataParser.processBarcodeMapping(thread);
        dataParser.writeLog(outFile);
        dataParser.writeAbCountsPerSc(outFile);

    }
    catch (const std::exception& e) 
    {
        std::cerr << "Counting failed with error: " << e.what() << "\n";
        return false;
    }

    return true;
}

